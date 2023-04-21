/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: kernels.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */


#ifndef _KERNELS_H_
#define _KERNELS_H_

#include <numeric>
#include <algorithm>
#include "miscellaneous.h"
#include "rotation.h"

#ifdef __CUDACC__
#include "helper_cuda.h"
#include <thrust/random.h>
#endif

#define GAMMA  267515315. // rad/s.T


#ifdef __CUDACC__
__host__  __device__ __forceinline__
#else
inline
#endif
uint64_t sub2ind(uint32_t x, uint32_t y, uint32_t z, uint32_t lenx, uint32_t leny)
{
    return (uint64_t(z*lenx*leny) + y*lenx + x);
}

template <typename T>
#ifdef __CUDACC__
__device__ __forceinline__
#endif
uint16_t find_ind(const T *arr, uint16_t len, T seek)
{
    for(uint16_t i=0; i<len; i++)
        if (arr[i] == seek)
            return i;
    return len;
}


#ifdef __CUDACC__
__global__ 
#endif
void simulation_kernel(const simulation_parameters *param, const float *pFieldMap, const bool *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint32_t spin_no = 0)
{
#ifdef __CUDACC__
    spin_no = blockIdx.x * blockDim.x + threadIdx.x ;
    if (spin_no >= param->n_spins)
        return;

    thrust::minstd_rand gen(param->seed + spin_no);
    thrust::normal_distribution<float> dist_random_walk_xyz(0.f, sqrt(6 * param->diffusion_const * param->dt));
    //thrust::uniform_real_distribution<float> dist_random_walk_xyz(-sqrt(6 * param.diffusion_const * param.dt), sqrt(6 * param.diffusion_const * param.dt));
    gen.discard(param->seed + spin_no);
#else
    std::minstd_rand  gen(param->seed + spin_no);
    std::normal_distribution<float> dist_random_walk_xyz(0.f, sqrt(6 * param->diffusion_const * param->dt));
#endif

    //uint16_t n_timepoints_local;
    float field = 0., phase_cycle = 0., time_elapsed = 0.; 
    float m0[3], m1[3]; 
    float xyz[3], xyz_new[3];
    for(uint32_t i=0, shift=3*spin_no; i<3; i++)
    {
        xyz[i] = XYZ0[shift + i];
        m0[i]  = M0[shift + i];
    }

    // -alpha/2 RF pulse (along x-axis) + TR/2 relaxation
    if (param->enApplyFA2)
    {
        //xrot(-param->s2, param->c2, m0, m1); // note this is -FA/2
        xrot_withphase (param->s2, param->c2, phase_cycle += param->phase_cycling, m0, m1);
        relax(param->e12, param->e22, m1, m0);
    }

    bool is_lastdummy = false;
    for (uint32_t dummy_scan = 0; dummy_scan < param->n_dummy_scan + 1; dummy_scan++)
    {
        is_lastdummy = (dummy_scan == param->n_dummy_scan);
        //n_timepoints_local = is_lastdummy ? (param->n_timepoints) / 2 : (param->n_timepoints); // random walk till TR/2 in the final execution of the loop
        
        //alpha RF pulse (along x-axis) + TR or TR/2 relaxation
        //float s_cycl = (dummy_scan % 2 == 0) ? param->s : -param->s; // PI phase cycling, starts with +FA (since we have -FA/2 above as the first pulse)
        //xrot(s_cycl, param->c, m0, m1);
        if (phase_cycle > 360.0)
            phase_cycle -= 360.0;
        else if (phase_cycle < 0)
            phase_cycle += 360.0;

        xrot_withphase (param->s, param->c, phase_cycle += param->phase_cycling, m0, m1);
        // copy m1 to m0
        for(uint8_t i=0; i<3; i++)
            m0[i] = m1[i];

        // random walk with boundries and accomulate phase
        uint64_t ind=0, ind_old=param->matrix_length+1;
        uint16_t current_timepoint = 0, old_timepoint = 0;
        float accumulated_phase = 0.f;

        while (current_timepoint < param->n_timepoints)
        {
            // random walk
            for (uint8_t i=0; i<3; i++)
            {
                xyz_new[i] = xyz[i] + dist_random_walk_xyz(gen); // new spin position after random-walk
                if (xyz_new[i] < 0)
                    xyz_new[i] += param->sample_length[i];
                else if (xyz_new[i] > param->sample_length[i])
                    xyz_new[i] -= param->sample_length[i];
            }
            // subscripts to linear indices
            ind = sub2ind(ROUND(xyz_new[0]*param->scale2grid[0]), ROUND(xyz_new[1]*param->scale2grid[1]), ROUND(xyz_new[2]*param->scale2grid[2]), param->fieldmap_size[0], param->fieldmap_size[1]);
            //accumulate phase
            if(ind != ind_old) // used this trick for fewer access to the global memory which is slow. Helpful for big samples!
            {               
                if (pMask[ind] == true) // check doesn't cross a vessel 
                    continue;       
                field = pFieldMap[ind_old = ind];
            }     
            accumulated_phase += field;

            // apply refocusing pulse if there is any
            if(param->enRefocusing)
            {
                uint16_t ind = find_ind<uint16_t>(param->T_SE, param->n_SE, current_timepoint);                
                if (ind < param->n_SE)
                {   
                    // dephase                 
                    accumulated_phase *= param->B0 * GAMMA * param->dt; // Fieldmap per Tesla to radian    
                    zrot(accumulated_phase * RAD2DEG, m0, m1); // dephase; zrot expects angle in degrees 
                    // relax
                    time_elapsed = (current_timepoint - old_timepoint) * param->dt;
                    relax(exp(-time_elapsed/param->T1), exp(-time_elapsed/param->T2), m1);
                    // refocusing pulse
                    xrot_withphase (param->RF_SE[ind], param->RF_SE_PHS[ind], m1, m0); // Note m0 and m1 are swapped here, so that we can use m0 for the next iteration
                    accumulated_phase = 0; // reset phase since we have applied it in the previous step
                    old_timepoint = current_timepoint;
                }
            }

            // record echo, only in the last scan
            if (is_lastdummy)
            {
                uint16_t ind = find_ind<uint16_t>(param->TE, param->n_TE, current_timepoint);
                if (ind < param->n_TE)
                {   
                    // dephase                  
                    accumulated_phase *= param->B0 * GAMMA * param->dt; // Fieldmap per Tesla to radian    
                    zrot(accumulated_phase * RAD2DEG, m0, m1); // dephase; zrot expects angle in degrees 
                    // relax                    
                    time_elapsed = (current_timepoint - old_timepoint) * param->dt;
                    relax(exp(-time_elapsed/param->T1), exp(-time_elapsed/param->T2), m1);
                    // save echo and copy m1 to m0 for the next iteration
                    for (uint32_t i=0, shift=3*param->n_TE*spin_no + 3*ind; i<3; i++)
                        M1[shift + i] = m0[i] = m1[i];   
                    accumulated_phase = 0; // reset phase since we have applied it in the previous step
                    old_timepoint = current_timepoint;
                }
            }

            for (uint8_t i = 0; i < 3; i++)
                xyz[i] = xyz_new[i];
 
            current_timepoint++;            
        }
        // dephase
        accumulated_phase *= param->B0 * GAMMA * param->dt; // Fieldmap per Tesla to radian     
        zrot(accumulated_phase * RAD2DEG, m0, m1); // dephase; zrot expects angle in degrees 
        // relax
        time_elapsed = (current_timepoint - old_timepoint) * param->dt;
        relax(exp(-time_elapsed/param->T1), exp(-time_elapsed/param->T2), m1);

        // copy m1 to m0 for the next iteration
        for(uint8_t i=0; i<3; i++)
            m0[i] = m1[i];
    }
    // save final position
    for (uint32_t i=0, shift=3*spin_no; i<3; i++)
        XYZ1[shift + i] = xyz[i];
}

#ifdef __CUDACC__
__global__ void scale_initial_positions(float *scaled_xyz, float *initial_xyz, float scale, uint32_t size)
{
    int n = blockIdx.x * blockDim.x + threadIdx.x ;
    if(n < size)
    {
        uint32_t ind = 3*n;
        scaled_xyz[ind+0] = initial_xyz[ind+0] * scale;
        scaled_xyz[ind+1] = initial_xyz[ind+1] * scale;
        scaled_xyz[ind+2] = initial_xyz[ind+2] * scale;
    }
}
#endif

// generate random initial position
#ifdef __CUDACC__
__global__ 
#endif
void generate_initial_position(float *spin_position_xyz, simulation_parameters *param, bool *pMask, uint32_t spin_no = 0)
{
#ifdef __CUDACC__
    spin_no = blockIdx.x * blockDim.x + threadIdx.x ;
    if(spin_no >= param->n_spins)
        return;

    thrust::minstd_rand  gen(param->seed + spin_no);
    thrust::uniform_real_distribution<float> dist_initial_point(0.1, 0.9);
    gen.discard(param->seed + spin_no);
#else
    std::minstd_rand  gen(param->seed + spin_no);
    std::uniform_real_distribution<float> dist_initial_point(0.1, 0.9);
 #endif

    float scale2grid[3];
    for(int i=0; i<3; i++)
        scale2grid[i] = (param->fieldmap_size[i]-1.) / param->sample_length[i];

    uint64_t index = 0;
    float *xyz = spin_position_xyz + 3*spin_no;
    do
    {
        for (uint8_t i = 0; i < 3; i++)
            xyz[i] = dist_initial_point(gen) * param->sample_length[i];
        index = sub2ind(ROUND(xyz[0]*scale2grid[0]), ROUND(xyz[1]*scale2grid[1]), ROUND(xyz[2]*scale2grid[2]), param->fieldmap_size[0], param->fieldmap_size[1]);
    }while (pMask[index] == true);
}


template <typename T>
uint32_t is_masked(std::vector<T> &XYZ0, std::vector<char> &mask, simulation_parameters *param)
{
    float scale2grid[3];
    for(uint8_t i=0; i<3; i++)
        scale2grid[i] = (param->fieldmap_size[i]-1.) / param->sample_length[i];

    std::vector<uint32_t> mask_counter(XYZ0.size()/3, false);
    #pragma omp parallel for
    for (int32_t i=0; i<XYZ0.size(); i+=3)
    {
        uint64_t index = sub2ind(ROUND(XYZ0[i] * scale2grid[0]), ROUND(XYZ0[i+1] * scale2grid[1]), ROUND(XYZ0[i+2] * scale2grid[2]), param->fieldmap_size[0], param->fieldmap_size[1]);
        mask_counter[i/3] = mask[index];
    }
    return std::accumulate(mask_counter.begin(), mask_counter.end(), 0);
}


void simulate_steady_state(simulation_parameters param)
{
    // alpha/2 RF pulse (along x-axis) + TE=TR/2 relaxation
    float m0_init[3] = { 0., 0., 1. };
    if (param.n_dummy_scan != 0)
    { // -alpha/2 RF pulse (along x-axis) + TE=TR/2 relaxation
        m0_init[1] = -sinf (-param.FA / 2.);
        m0_init[2] = cosf (-param.FA / 2.);
        relax (param.e12, param.e22, m0_init);
    }

    float m0t[3], m1t[3];
    std::copy (m0_init, m0_init + 3, m0t);
    std::cout << "Steady-state report:" << std::endl;
    std::cout << "\tM0 after alpha/2 pulse & TR/2 relaxation = [" << m0t[0]
              << " " << m0t[1] << " " << m0t[2] << "]" << std::endl;
    float phase_cycle = 0;
    for (uint32_t rep = 0; rep < 3; rep++)
        for (uint32_t dummy_scan = 0; dummy_scan < param.n_dummy_scan; dummy_scan++)
        {
            // s_cycl = (dummy_scan % 2 == 0) ? -param.s : param.s;
            if (phase_cycle > 360.0)
                phase_cycle -= 360.0;
            else if (phase_cycle < 0)
                phase_cycle += 360.0;
            xrot_withphase (param.s, param.c, phase_cycle += param.phase_cycling, m0t, m1t);

            if (dummy_scan != param.n_dummy_scan - 1)
                relax (param.e1, param.e2, m1t);
            else
            {
                relax (param.e12, param.e22, m1t);
                std::cout << "\tMagnetization after "
                          << param.n_dummy_scan * (rep + 1) << " RF shots = ["
                          << m1t[0] << " " << m1t[1] << " " << m1t[2] << "]"
                          << std::endl;
                relax (param.e12, param.e22, m1t);
            }
            std::copy (m1t, m1t + 3, m0t);
        }
}

#endif // _KERNELS_H_