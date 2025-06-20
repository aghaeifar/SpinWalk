/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: kernels.cu
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : 
 * -------------------------------------------------------------------------- */

#include <cinttypes>
#include <algorithm>
#include <execution>
#include "kernels.cuh"
#include "definitions.h"

#ifdef __CUDACC__
#include "helper_cuda.h"
#include <cuda_runtime.h>
#include <thrust/random.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#endif

#define ABS(x) ((x) < 0 ? -(x) : (x))

namespace sim {

uint8_t find_max(const std::vector<uint8_t> &data)
{
#ifdef __CUDACC__
    thrust::device_vector<uint8_t> thrust_vec(data.begin(), data.end());
    uint8_t m = *thrust::max_element(thrust_vec.begin(), thrust_vec.end());
#else
    uint8_t m = *std::max_element(std::execution::par, data.begin(), data.end());
#endif
    return m;
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void dephase_relax(float *m0, float *m1, float accumulated_phase, float T1, float T2, float time_elapsed)
{
    // dephase                
    zrot(accumulated_phase, m0, m1); 
    // relax
    if (T1 >= 0 && T2 >= 0)   
        relax(exp(-time_elapsed/T1), exp(-time_elapsed/T2), m1);
}


#ifdef __CUDACC__
__global__ 
void cu_sim(const parameters param, const parameters_uvec param_uvec, const float *pFieldMap, const uint8_t *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint8_t *T)
{
    uint32_t spin_no = blockIdx.x * blockDim.x + threadIdx.x ;
    if (spin_no >= param.n_spins)
        return;  
    sim(param, param_uvec, pFieldMap,pMask,M0, XYZ0, M1, XYZ1, T, spin_no);
}
#endif 




#ifdef __CUDACC__
__host__  __device__ 
#endif 
void sim(const parameters &param, const parameters_uvec &param_uvec, const float *pFieldMap, const uint8_t *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint8_t *T, uint32_t spin_no)
{    
    // printf("spin=%d\n", spin_no);
    float *xyz1 = XYZ1 + 3*spin_no * (param.enRecordTrajectory ? (param.n_dummy_scan + 1)*(param.n_timepoints) : 1);
#ifdef __CUDACC__
    thrust::minstd_rand gen_r(param.seed + spin_no);
    thrust::minstd_rand gen_u(param.seed + spin_no);
    thrust::normal_distribution<float> dist_random_walk_xyz(0.f, 1.0f);
    thrust::uniform_real_distribution<float> dist_cross_tissue(0.f, 1.f);
#else
    std::mt19937 gen_r(param.seed + spin_no);
    std::mt19937 gen_u(param.seed + spin_no);
    std::normal_distribution<float> dist_random_walk_xyz(0.f, 1.0f);
    std::uniform_real_distribution<float> dist_cross_tissue(0.f, 1.f);
#endif
    gen_r.discard(param.seed + spin_no); // each spins has its own seed, and param.seed differes for each GPU in HPC with multiple GPUs
    gen_u.discard(param.seed + spin_no); // each spins has its own seed, and param.seed differes for each GPU in HPC with multiple GPUs
    
    uint32_t itr = 0;
    float field = 0., T1=0., T2=0., rf_phase = param_uvec.RF_PH_deg.ptr[0], time_elapsed = 0.; 
    float m0[3], m1[3]; 
    double xyz_old[3], xyz_new[3], scale2grid[3];

    for(uint32_t i=0, shift=3*spin_no; i<3; i++) {
        xyz_old[i] = xyz_new[i] = xyz1[i] = XYZ0[shift + i];
        m0[i] = M0[shift + i];
        scale2grid[i] = param.phantom_size[i] / param.fov[i];
    }
    // tissue type
    uint8_t ts, ts_old;    
    auto indx = sub2ind(xyz1[0]*scale2grid[0], xyz1[1]*scale2grid[1], xyz1[2]*scale2grid[2], param.phantom_size[0], param.phantom_size[1], param.phantom_size[2]);
    ts = ts_old = pMask[indx];
    double diffusivity_scale = param_uvec.diffusivity.ptr[ts_old];
    
    bool is_lastscan = false;
    for (uint32_t dummy_scan = 0; dummy_scan < param.n_dummy_scan + 1; dummy_scan++) {
        is_lastscan = (dummy_scan == param.n_dummy_scan);
        
        float new_rf_phase = rf_phase + dummy_scan*param.linear_phase_cycling + dummy_scan*(dummy_scan+1)/2.0*param.quadratic_phase_cycling;
        while (new_rf_phase > 360.0)
            new_rf_phase -= 360.0;
        while (new_rf_phase < 0)
            new_rf_phase += 360.0;
        
        // ------ apply the first RF pulse. The start time for the first RF pulse is always 0 ------
        xrot_withphase (param.s, param.c, new_rf_phase, m0, m1);

        for(uint8_t i=0; i<3; i++) // copy m1 to m0
            m0[i] = m1[i];

        // define and reset all counters
        int64_t ind=0, ind_old=param.matrix_length+1;
        uint32_t current_timepoint = 0, old_timepoint = 0;
        uint16_t current_rf = 1, current_te = 0, counter_dephasing = 0, counter_gradient = 0;
        float accumulated_phase = 0.f;      
        // ------ loop over timepoints ------
        while (current_timepoint < param.n_timepoints) { // param.n_timepoints is the total number of timepoints (= TR/dwelltime)
            // ------ generate random walks and wrap around the boundries ------
            for (uint8_t i=0; i<3 && diffusivity_scale != 0.; i++) {
                double rnd_wlk = dist_random_walk_xyz(gen_r) * diffusivity_scale;
                xyz_new[i] = xyz_old[i] + rnd_wlk; // new spin position after random-walk
                if (xyz_new[i] < 0)
                    xyz_new[i] += (param.enCrossFOV ? param.fov[i] : 2*std::abs(rnd_wlk)); // rnd_wlk is negative here
                else if (xyz_new[i] >= param.fov[i])
                    xyz_new[i] -= (param.enCrossFOV ? param.fov[i] : 2*std::abs(rnd_wlk)); // rnd_wlk is positive here
            }
            
            // ------ subscripts to linear indices ------
            ind = sub2ind(xyz_new[0]*scale2grid[0], xyz_new[1]*scale2grid[1], xyz_new[2]*scale2grid[2], param.phantom_size[0], param.phantom_size[1], param.phantom_size[2]);
            if(ind >= param.matrix_length || ind < 0) {
                printf("\n--------------------- <Error> ---------------------\n");
                printf("spin = %d\ntimepoint = %d\nind = %" PRId64 "\nMatrixSize = %" PRId64 "\nPhantomSize = (%" PRId64 ", %" PRId64 ", %" PRId64 ")\n", spin_no, current_timepoint, ind, param.matrix_length, param.phantom_size[0], param.phantom_size[1], param.phantom_size[2]);
                printf("FoV     = (%.10f, %.10f, %.10f)\nxyz_new = (%.10f, %.10f, %.10f)\nxyz_old = (%.10f, %.10f, %.10f)\n", param.fov[0], param.fov[1], param.fov[2], xyz_new[0], xyz_new[1], xyz_new[2], xyz_old[0], xyz_old[1], xyz_old[2]);
                printf("Error = (%d, %d, %d)\nscale2grid = (%.5f, %.5f, %.5f)\nscale2grid*xyz_new = (%.10f, %.10f, %.10f)\n", xyz_new[0] >= param.fov[0], xyz_new[1] >= param.fov[1], xyz_new[2] >= param.fov[2], scale2grid[0], scale2grid[1], scale2grid[2], scale2grid[0]*xyz_new[0], scale2grid[1]*xyz_new[1], scale2grid[2]*xyz_new[2]);
                return;
            }
              
            // ------ accumulate phase ------
            if(ind != ind_old) {  // fewer access to the global memory which is slow. Helpful for large samples!
                // cross-tissue diffusion
                ts = pMask[ind];          
                if (ts != ts_old) {
                    if (dist_cross_tissue(gen_u) >= param_uvec.pXY.ptr[ts_old*param.n_substrate + ts]) {
                        if(itr++ > param.max_iterations) {
                            printf("Warning! spin %d is stuck at (%f, %f, %f) and is considered lost (dummy=%d time=%d).\n", spin_no, xyz_new[0], xyz_new[1], xyz_new[2], dummy_scan, current_timepoint);
                            // we must reset magnetization for all echoes here! 
                            return;
                        }
                        continue;
                    }
                    ts_old = ts;
                }
                    
                ind_old = ind;   
                field = pFieldMap != nullptr ? pFieldMap[ind]:0.f;
                T1 = param_uvec.T1_ms.ptr[ts_old] * 1e-3; // ms -> s
                T2 = param_uvec.T2_ms.ptr[ts_old] * 1e-3; // ms -> s
                diffusivity_scale = param_uvec.diffusivity.ptr[ts_old];
            }   
            accumulated_phase += field;
            itr = 0;

            // ------ apply ideal dephasing if there is any ------
            if(counter_dephasing < param_uvec.dephasing_us.size && param_uvec.dephasing_us.ptr[counter_dephasing] == current_timepoint) {
                accumulated_phase += (float)spin_no * param_uvec.dephasing_deg.ptr[counter_dephasing] / (float)param.n_spins; // assign dephasing linearly to spins 
                counter_dephasing++;
            }

            // ------ apply gradient if there is any ------
            if(counter_gradient < param_uvec.gradient_us.size && param_uvec.gradient_us.ptr[counter_gradient] == current_timepoint) {
                const float Gx = *(param_uvec.gradientX_mTm.ptr + counter_gradient);
                const float Gy = *(param_uvec.gradientY_mTm.ptr + counter_gradient);
                const float Gz = *(param_uvec.gradientZ_mTm.ptr + counter_gradient);
                accumulated_phase += (Gx*xyz_new[0] + Gy*xyz_new[1] + Gz*xyz_new[2]) * 1e-3 * param.timestep_us * 1e-6 * GAMMA * RAD2DEG; //  Gx * x + Gy * y + Gz * z
                counter_gradient++;
            }
                
            // ------ apply other RF pulse if there is any ------
            if(current_rf < param_uvec.RF_us.size && param_uvec.RF_us.ptr[current_rf] == current_timepoint) {
                // dephase and relax    
                time_elapsed = (current_timepoint - old_timepoint) * param.timestep_us * 1e-6;
                dephase_relax(m0, m1, accumulated_phase, T1, T2, time_elapsed);
                // apply RF pulse
                xrot_withphase (param_uvec.RF_FA_deg.ptr[current_rf], param_uvec.RF_PH_deg.ptr[current_rf], m1, m0); // Note m0 and m1 are swapped here, so that we can use m0 for the next iteration
                accumulated_phase = 0; // reset phase since we have it now applied
                old_timepoint = current_timepoint;
                current_rf++;
            }

            // ------ echoes are only recorded in the last scan ------
            if (is_lastscan && current_te < param_uvec.TE_us.size && param_uvec.TE_us.ptr[current_te] == current_timepoint) {
                // dephase and relax                
                time_elapsed = (current_timepoint - old_timepoint) * param.timestep_us * 1e-6;
                dephase_relax(m0, m1, accumulated_phase, T1, T2, time_elapsed);
                // save echo and copy m1 to m0 for the next iteration
                for (uint32_t i=0, shift=3*param_uvec.TE_us.size*spin_no + 3*current_te; i<3; i++)
                    M1[shift + i] = m0[i] = m1[i];
                
                T[spin_no*param_uvec.TE_us.size + current_te] = ts_old;

                accumulated_phase = 0; // reset phase since we have applied it in the previous step
                old_timepoint = current_timepoint;
                current_te++;
            }

            // update old position with the new one
            if(param.enRecordTrajectory && (current_timepoint != 0 || dummy_scan != 0))
                xyz1 += 3;
            for (uint8_t i=0; i < 3; i++)
                xyz1[i] = xyz_old[i] = xyz_new[i];
            // increase timepoint
            current_timepoint++;            
        }
        // dephase and relax    
        time_elapsed = (current_timepoint - old_timepoint) * param.timestep_us * 1e-6;
        dephase_relax(m0, m1, accumulated_phase, T1, T2, time_elapsed);

        // copy m1 to m0 for the next iteration
        for(uint8_t i=0; i<3; i++)
            m0[i] = m1[i];
    }
}

} // namespace sim