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
int64_t sub2ind(int32_t x, int32_t y, int32_t z, int32_t lenx, int32_t leny)
{
    return (int64_t(z*lenx*leny) + y*lenx + x);
}

#ifdef __CUDACC__
__global__ 
#endif
void simulation_kernel(simulation_parameters *param, float *pFieldMap, bool *pMask, float *random_walk_xyz_init_scaled, float *mxyz, int32_t spin_no = -1)
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

    int16_t n_timepoints_local;
    float accumulated_phase = 0., field = 0.; 
    float m0[3] = {0., 0., 1.}, m1[3]; 
    float xyz[3], xyz_new[3];
    for(int32_t i=0; i<3; i++)
        xyz[i] = random_walk_xyz_init_scaled[3*spin_no + i];

    // -alpha/2 RF pulse (along x-axis) + TR/2 relaxation
    if (param->n_dummy_scan != 0)
    {
        m0[1] = param->s2; // note this is -FA/2 : +param->s2
        m0[2] = param->c2;
        relax(param->e12, param->e22, m0);
    }

    bool is_lastdummy = false;
    for (int32_t dummy_scan = 0; dummy_scan < param->n_dummy_scan + 1; dummy_scan++)
    {
        is_lastdummy = (dummy_scan == param->n_dummy_scan);
        n_timepoints_local = is_lastdummy ? (param->n_timepoints) / 2 : (param->n_timepoints); // random walk till TR/2 in the final execution of the loop
        
        // alpha RF pulse (along x-axis) + TR or TR/2 relaxation
        float s_cycl = (dummy_scan % 2 == 0) ? param->s : -param->s; // PI phase cycling, starts with +FA (since we have -FA/2 above as the first pulse)
        xrot(s_cycl, param->c, m0, m1);

        // copy m1 to m0
        for(int32_t i=0; i<3; i++)
            m0[i] = m1[i];

        // random walk with boundries and accomulate phase
        int64_t ind=0, ind_old=-1;
        int16_t current_timepoint = 0;
        accumulated_phase = 0;

        while (current_timepoint < n_timepoints_local)
        {
            for (int32_t i=0; i<3; i++)
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
            if(ind != ind_old) // used this trick to access fewer to the global memorsy which is slow. Helpful for big samples!
            {               
                if (pMask[ind] == true) // check doesn't cross a vessel 
                    continue;       
                field = pFieldMap[ind];
                ind_old = ind; 
            }
     
            accumulated_phase += field;
            if(param->enRefocusing180 && current_timepoint == param->n_timepoints/4)
                accumulated_phase *= -1; // 180 degree refocusing pulse

            for (int32_t i = 0; i < 3; i++)
                xyz[i] = xyz_new[i];
 
            current_timepoint++;            
        }
        
        accumulated_phase *= param->B0 * GAMMA * param->dt; // Fieldmap per Tesla to radian     
        zrot(accumulated_phase, m0, m1); // dephase
        relax(is_lastdummy ? param->e12 : param->e1, is_lastdummy ? param->e22 : param->e2, m1);

        // copy m1 to m0 for the next iteration
        for(int32_t i=0; i<3; i++)
            m0[i] = m1[i];
    }

    int32_t ind = 3*spin_no;
    for (int32_t i=0; i<3; i++)
        mxyz[ind + i] = m1[i];

}

#ifdef __CUDACC__
__global__ void scale_initial_positions(float *scaled_xyz, float *initial_xyz, float scale, int32_t size)
{
    int n = blockIdx.x * blockDim.x + threadIdx.x ;
    if(n < size)
    {
        int32_t ind = 3*n;
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
void generate_initial_position(float *spin_position_xyz, simulation_parameters *param, bool *pMask, int32_t spin_no = -1)
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

    int32_t index = 0;
    float *xyz = spin_position_xyz + 3*spin_no;
    do
    {
        for (int32_t i = 0; i < 3; i++)
            xyz[i] = dist_initial_point(gen) * param->sample_length[i];
        index = sub2ind(ROUND(xyz[0]*scale2grid[0]), ROUND(xyz[1]*scale2grid[1]), ROUND(xyz[2]*scale2grid[2]), param->fieldmap_size[0], param->fieldmap_size[1]);
    }while (pMask[index] == true);
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

    float m0t[3], m1t[3], s_cycl;
    std::copy (m0_init, m0_init + 3, m0t);
    std::cout << "Steady-state report:" << std::endl;
    std::cout << "\tM0 after alpha/2 pulse & TR/2 relaxation = [" << m0t[0]
              << " " << m0t[1] << " " << m0t[2] << "]" << std::endl;
    for (int32_t rep = 0; rep < 3; rep++)
        for (int32_t dummy_scan = 0; dummy_scan < param.n_dummy_scan; dummy_scan++)
        {
            s_cycl = (dummy_scan % 2 == 0) ? -param.s : param.s;
            xrot (s_cycl, param.c, m0t, m1t);
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