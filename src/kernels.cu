/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: kernels.cuh
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */


#include <algorithm>
#include "kernels.cuh"
#include "rotation.cuh"
#include "helper_cuda.h"
#include <cuda_runtime.h>
#include <boost/log/trivial.hpp>
#include <thrust/random.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------

uint8_t find_max(const std::vector<uint8_t> &data)
{
    thrust::device_vector<uint8_t> gpu_vec(data.begin(), data.end());
    uint8_t m = *thrust::max_element(gpu_vec.begin(), gpu_vec.end());
    return m;
}

__device__ __forceinline__ void dephase_relax(float *m0, float *m1, float accumulated_phase, float T1, float T2, float time_elapsed)
{
    // dephase                
    zrot(accumulated_phase, m0, m1); 
    // relax
    relax(exp(-time_elapsed/T1), exp(-time_elapsed/T2), m1);
}

__global__ void cu_sim(const simulation_parameters *param, const float *pFieldMap, const uint8_t *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint8_t *T)
{
    auto spin_no = blockIdx.x * blockDim.x + threadIdx.x ;
    if (spin_no >= param->n_spins)
        return;

    float *xyz1 = XYZ1 + 3*spin_no * (param->enRecordTrajectory ? (param->n_dummy_scan + 1)*(param->n_timepoints) : 1);

    thrust::minstd_rand gen_r(param->seed + spin_no);
    thrust::minstd_rand gen_u(param->seed + spin_no);
    thrust::normal_distribution<float> dist_random_walk_xyz(0., sqrt(6. * param->diffusion_const * param->dt));
    // thrust::uniform_real_distribution<float> dist_random_walk_xyz(-sqrt(6 * param->diffusion_const * param->dt), sqrt(6 * param->diffusion_const * param->dt));
    thrust::uniform_real_distribution<float> dist_cross_tissue(0.0f, 1.0f);
    gen_r.discard(param->seed + spin_no); // each spins has its own seed, and param->seed differes for each GPU in HPC with multiple GPUs
    gen_u.discard(param->seed + spin_no); // each spins has its own seed, and param->seed differes for each GPU in HPC with multiple GPUs

    uint32_t itr = 0;
    float field = 0., T1=0., T2=0., rf_phase = param->RF_PH[0], time_elapsed = 0.; 
    float m0[3], m1[3]; 
    float xyz_new[3];
    for(uint32_t i=0, shift=3*spin_no; i<3; i++)
    {
        xyz1[i] = XYZ0[shift + i];
        m0[i]  = M0[shift + i];
    }
    // tissue type
    uint8_t ts, ts_old;
    auto indx = sub2ind(xyz1[0]*param->scale2grid[0], xyz1[1]*param->scale2grid[1], xyz1[2]*param->scale2grid[2], param->fieldmap_size[0], param->fieldmap_size[1]);
    ts_old = pMask[indx];

    bool is_lastscan = false;
    for (uint32_t dummy_scan = 0; dummy_scan < param->n_dummy_scan + 1; dummy_scan++)
    {
        is_lastscan = (dummy_scan == param->n_dummy_scan);
        
        while (rf_phase > 360.0)
            rf_phase -= 360.0;
        while (rf_phase < 0)
            rf_phase += 360.0;
        
        // ------ apply the first RF pulse. The start time for the first RF pulse is always 0 ------
        xrot_withphase (param->s, param->c, rf_phase += param->phase_cycling, m0, m1);

        for(uint8_t i=0; i<3; i++) // copy m1 to m0
            m0[i] = m1[i];

        // ------ loop over timepoints ------
        uint64_t ind=0, ind_old=param->matrix_length+1;
        uint32_t current_timepoint = 0, old_timepoint = 0;
        uint16_t current_rf = 1, current_te = 0, counter_dephasing = 0, counter_gradient = 0;
        float accumulated_phase = 0.f;        
        while (current_timepoint < param->n_timepoints) // param->n_timepoints is the total number of timepoints (= TR/dwelltime)
        {
            // ------ generate random walks and wrap around the boundries ------
            float rnd_wlk;
            for (uint8_t i=0; i<3; i++)
            {
                rnd_wlk = dist_random_walk_xyz(gen_r);
                xyz_new[i] = xyz1[i] + rnd_wlk; // new spin position after random-walk
                if (xyz_new[i] < 0)
                    xyz_new[i] += param->enCrossBoundry ? param->sample_length[i] : -2*rnd_wlk; // rnd_wlk is negative here
                else if (xyz_new[i] > param->sample_length[i])
                    xyz_new[i] -= param->enCrossBoundry ? param->sample_length[i] : 2*rnd_wlk;
            }
           
            // ------ subscripts to linear indices ------
            ind = sub2ind(xyz_new[0]*param->scale2grid[0], xyz_new[1]*param->scale2grid[1], xyz_new[2]*param->scale2grid[2], param->fieldmap_size[0], param->fieldmap_size[1]);
            if(ind > param->matrix_length)
            {
                printf("Error:spin=%d, ind=%llu, %d,  scale=(%f), xyz_new=(%f, %f, %f)\n",spin_no, ind, current_timepoint, param->scale2grid[0], xyz_new[0], xyz_new[1], xyz_new[2]);
                return;
            }
            // ------ accumulate phase ------
            if(ind != ind_old) // fewer access to the global memory which is slow. Helpful for large samples!
            {   
                // cross-tissue diffusion
                ts = pMask[ind];           
                if (ts != ts_old) 
                    if (dist_cross_tissue(gen_u) >= param->pXY[ts_old*param->n_tissue_type + ts])
                    {
                        if(itr++ > param->max_iterations)
                        {
                            printf("Warning! spin %d is stuck at (%f, %f, %f) and is considered lost (dummy=%d time=%d).\n", spin_no, xyz_new[0], xyz_new[1], xyz_new[2], dummy_scan, current_timepoint);
                            return;
                        }
                        continue;
                    }
                    else
                        ts_old = ts;
                itr = 0;       
                field = pFieldMap[ind_old = ind];
                ind = pMask[ind]; // the index of the tissue type
                T1 = param->T1[ind];
                T2 = param->T2[ind];
            }     
            accumulated_phase += field;
         
            // ------ apply dephasing if there is any ------
            if(counter_dephasing < param->n_dephasing && param->dephasing_T[counter_dephasing] == current_timepoint)
            {
                accumulated_phase += (float)spin_no * param->dephasing[counter_dephasing] / (float)param->n_spins; // assign dephasing linearly to spins 
                counter_dephasing++;
            }

            // ------ apply gradient if there is any ------
            if(counter_gradient < param->n_gradient && param->gradient_T[counter_gradient] == current_timepoint)
            {
                const float *Gxyz = param->gradient_xyz + 3*counter_gradient;
                accumulated_phase +=  (Gxyz[0]*xyz_new[0] + Gxyz[1]*xyz_new[1] + Gxyz[2]*xyz_new[2]) * param->dt*GAMMA*RAD2DEG; //  Gx * x + Gy * y + Gz * z
                counter_gradient++;
            }
                 
            // ------ apply other RF pulse if there is any ------
            if(current_rf < param->n_RF && param->RF_ST[current_rf] == current_timepoint)
            {
                // dephase and relax    
                time_elapsed = (current_timepoint - old_timepoint) * param->dt;
                dephase_relax(m0, m1, accumulated_phase, T1, T2, time_elapsed);
                // apply RF pulse
                xrot_withphase (param->RF_FA[current_rf], param->RF_PH[current_rf], m1, m0); // Note m0 and m1 are swapped here, so that we can use m0 for the next iteration
                accumulated_phase = 0; // reset phase since we have it now applied
                old_timepoint = current_timepoint;
                current_rf++;
            }

            // ------ echoes are only recorded in the last scan ------
            if (is_lastscan && current_te < param->n_TE && param->TE[current_te] == current_timepoint)
            {
                // dephase and relax                
                time_elapsed = (current_timepoint - old_timepoint) * param->dt;
                dephase_relax(m0, m1, accumulated_phase, T1, T2, time_elapsed);
                // save echo and copy m1 to m0 for the next iteration
                for (uint32_t i=0, shift=3*param->n_TE*spin_no + 3*current_te; i<3; i++)
                    M1[shift + i] = m0[i] = m1[i];
                T[spin_no*param->n_TE + current_te] = ts;

                accumulated_phase = 0; // reset phase since we have applied it in the previous step
                old_timepoint = current_timepoint;
                current_te++;
            }

            // update old position with the new one
            if(param->enRecordTrajectory && (current_timepoint != 0 || dummy_scan != 0))
                xyz1 += 3;      
            for (uint8_t i=0; i < 3; i++)
                xyz1[i] = xyz_new[i];
            // increase timepoint
            current_timepoint++;            
        }
        // dephase and relax    
        time_elapsed = (current_timepoint - old_timepoint) * param->dt;
        dephase_relax(m0, m1, accumulated_phase, T1, T2, time_elapsed);

        // copy m1 to m0 for the next iteration
        for(uint8_t i=0; i<3; i++)
            m0[i] = m1[i];
    }
}


//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------

__global__ void cu_scalePos(float *scaled_xyz, float *initial_xyz, float scale, uint64_t size)
{
    uint64_t n = blockIdx.x * blockDim.x + threadIdx.x ;
    if(n < size)
    {
        uint64_t ind = 3*n;
        scaled_xyz[ind+0] = initial_xyz[ind+0] * scale;
        scaled_xyz[ind+1] = initial_xyz[ind+1] * scale;
        scaled_xyz[ind+2] = initial_xyz[ind+2] * scale;
    }
}

//---------------------------------------------------------------------------------------------
// CUDA kernel to perform array multiplication with a constant
//---------------------------------------------------------------------------------------------
__global__ void cu_scaleArray(float *array, float scale, uint64_t size)
{
    auto n = blockIdx.x * blockDim.x + threadIdx.x ;
    if(n < size)
        array[n] *= scale;
}

//---------------------------------------------------------------------------------------------
// CUDA kernel to generate random initial position
//---------------------------------------------------------------------------------------------

__global__ void cu_randPosGen(float *spin_position_xyz, simulation_parameters *param, const uint8_t *pMask, uint32_t spin_no)
{
    spin_no = blockIdx.x * blockDim.x + threadIdx.x ;
    if(spin_no >= param->n_spins)
        return;

    thrust::minstd_rand  gen(param->seed + spin_no);
    thrust::uniform_real_distribution<float> dist_initial_point(0.f, 1.f);
    gen.discard(param->seed + spin_no);

    float *xyz = spin_position_xyz + 3*spin_no;
    for (uint8_t i = 0; i < 3; i++)
        xyz[i] = dist_initial_point(gen) * param->sample_length[i];
}

//---------------------------------------------------------------------------------------------
//  check for CUDA and GPU device
//---------------------------------------------------------------------------------------------
void print_device_info()
{
    const int kb = 1024;
    const int mb = kb * kb;
    size_t free, total;
    
    int32_t cuda_version, driver_version;
    int32_t device_count = getDeviceCount();
    cudaRuntimeGetVersion(&cuda_version);
    cudaDriverGetVersion(&driver_version);
    std::cout << "\nDriver version: "<< driver_version << ", CUDA version: "<< cuda_version << ", Number of devices: " << device_count << std::endl;

    cudaDeviceProp device_properties;
    for(int i=0; i<device_count; i++)
    {
        cudaSetDevice(i);
        cudaGetDeviceProperties(&device_properties, i);
        std::cout << "Device " << i+1 << ", " << device_properties.name << std::endl;
        // cudaDeviceReset();
        cudaMemGetInfo(&free, &total);
        std::cout << "Free GPU memory: " << free / mb << " MB (out of " << total / mb << " MB)" << std::endl;
    }
}

uint32_t getDeviceCount()
{
    int32_t device_count;
    checkCudaErrors(cudaGetDeviceCount(&device_count));
    return device_count;
}

bool check_memory_size(size_t required_size_MB)
{
    size_t free, total;
    bool memory_ok = true;
    int32_t device_count = getDeviceCount();
    cudaDeviceProp device_properties;
    std::cout << "Device(s) memeory check:" << '\n';
    for(int i=0; i<device_count; i++)
    {
        cudaSetDevice(i);
        cudaGetDeviceProperties(&device_properties, i);
        cudaMemGetInfo(&free, &total);
        std::cout << "  Device " << i+1 << ", " << device_properties.name  << ": " << (free>required_size_MB ? "OK" : "Not enough") << '\n';
        if(free<required_size_MB)
            BOOST_LOG_TRIVIAL(fatal) << "Not enough GPU memory for the simulation in device "<< i <<"! Required=" << required_size_MB <<" MB, Available=" << free << " MB";
        memory_ok = free<required_size_MB ? false:memory_ok;
    }
    return memory_ok;
}
