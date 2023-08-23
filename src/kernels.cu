/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: kernels.cuh
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */


#include <algorithm>
#include <thrust/random.h>
#include "kernels.cuh"
#include "rotation.cuh"
#include "helper_cuda.h"

#define GAMMA  267515315. // rad/s.T


//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------

__global__ void cu_sim(const simulation_parameters *param, const float *pFieldMap, const bool *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint32_t spin_no)
{
    spin_no = blockIdx.x * blockDim.x + threadIdx.x ;
    if (spin_no >= param->n_spins)
        return;

    thrust::minstd_rand gen(param->seed + spin_no);
    thrust::normal_distribution<float> dist_random_walk_xyz(0.f, sqrt(6 * param->diffusion_const * param->dt));
    //thrust::uniform_real_distribution<float> dist_random_walk_xyz(-sqrt(6 * param.diffusion_const * param.dt), sqrt(6 * param.diffusion_const * param.dt));
    gen.discard(param->seed + spin_no); // each spins has its own seed, but this may not be true in an HPC with multiple GPUs

    //uint16_t n_timepoints_local;
    float field = 0., rf_phase = 0., time_elapsed = 0.; 
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
        xrot_withphase (param->s2, param->c2, rf_phase += param->phase_cycling, m0, m1);
        relax(param->e12, param->e22, m1, m0);
    }

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
        uint16_t current_timepoint = 0, old_timepoint = 0, current_rf = 1, current_te = 0, counter_dephasing = 0;
        float accumulated_phase = 0.f, dephase_deg = 0.f;
        while (current_timepoint < param->n_timepoints) // param->n_timepoints is the total number of timepoints (= TR/dwelltime)
        {
            // ------ generate random walks and wrap around the boundries ------
            for (uint8_t i=0; i<3; i++)
            {
                xyz_new[i] = xyz[i] + dist_random_walk_xyz(gen); // new spin position after random-walk
                if (xyz_new[i] < 0)
                    xyz_new[i] += param->sample_length[i];
                else if (xyz_new[i] > param->sample_length[i])
                    xyz_new[i] -= param->sample_length[i];
            }

            // ------ subscripts to linear indices ------
            ind = sub2ind(ROUND(xyz_new[0]*param->scale2grid[0]), ROUND(xyz_new[1]*param->scale2grid[1]), ROUND(xyz_new[2]*param->scale2grid[2]), param->fieldmap_size[0], param->fieldmap_size[1]);
            
            // ------ accumulate phase ------
            if(ind != ind_old) // used this trick for fewer access to the global memory which is slow. Helpful for large samples!
            {               
                if (pMask[ind] == true) // check doesn't cross a vessel 
                    continue;       
                field = pFieldMap[ind_old = ind];
            }     
            accumulated_phase += field;

            // ------ apply dephasing if there is any ------
            if(counter_dephasing < param->n_dephasing && param->dephasing_T[counter_dephasing] == current_timepoint)
            {
                float dephase_deg  =  (float)spin_no * param->dephasing[counter_dephasing] / (float)param->n_spins; // assign dephasing linearly to spins 
                accumulated_phase +=  dephase_deg / (param->B0 * GAMMA * param->dt * RAD2DEG); // scale dephasing to Tesla per dt
                counter_dephasing++;
            }
                 

            // ------ apply other RF pulse if there is any ------
            if(current_rf < param->n_RF && param->RF_ST[current_rf] == current_timepoint)
            {
                // dephase                
                dephase_deg =  accumulated_phase * param->B0 * GAMMA * param->dt * RAD2DEG; // convert accumulated phase to degree
                zrot(dephase_deg, m0, m1); 
                // relax
                time_elapsed = (current_timepoint - old_timepoint) * param->dt;
                relax(exp(-time_elapsed/param->T1), exp(-time_elapsed/param->T2), m1);
                // apply RF pulse
                xrot_withphase (param->RF_FA[current_rf], param->RF_PH[current_rf], m1, m0); // Note m0 and m1 are swapped here, so that we can use m0 for the next iteration
                accumulated_phase = 0; // reset phase since we have it now applied
                old_timepoint = current_timepoint;
                current_rf++;
            }

            // ------ echoes are only recorded in the last scan ------
            if (is_lastscan && current_te < param->n_TE && param->TE[current_te] == current_timepoint)
            {
                // dephase                  
                dephase_deg =  accumulated_phase * param->B0 * GAMMA * param->dt * RAD2DEG; // convert accumulated phase to degree  
                zrot(dephase_deg, m0, m1);
                // relax                    
                time_elapsed = (current_timepoint - old_timepoint) * param->dt;
                relax(exp(-time_elapsed/param->T1), exp(-time_elapsed/param->T2), m1);
                // save echo and copy m1 to m0 for the next iteration
                for (uint32_t i=0, shift=3*param->n_TE*spin_no + 3*current_te; i<3; i++)
                    M1[shift + i] = m0[i] = m1[i];   
                accumulated_phase = 0; // reset phase since we have applied it in the previous step
                old_timepoint = current_timepoint;
                current_te++;
            }

            // update old position with the new one
            for (uint8_t i = 0; i < 3; i++)
                xyz[i] = xyz_new[i];
 
            current_timepoint++;            
        }
        // dephase
        dephase_deg =  accumulated_phase * param->B0 * GAMMA * param->dt * RAD2DEG; // convert accumulated phase to degree     
        zrot(dephase_deg, m0, m1);
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


//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------

__global__ void cu_scalePos(float *scaled_xyz, float *initial_xyz, float scale, uint32_t size)
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


//---------------------------------------------------------------------------------------------
// generate random initial position
//---------------------------------------------------------------------------------------------

__global__ void cu_randPosGen(float *spin_position_xyz, simulation_parameters *param, bool *pMask, uint32_t spin_no)
{
    spin_no = blockIdx.x * blockDim.x + threadIdx.x ;
    if(spin_no >= param->n_spins)
        return;

    thrust::minstd_rand  gen(param->seed + spin_no);
    thrust::uniform_real_distribution<float> dist_initial_point(0.1, 0.9);
    gen.discard(param->seed + spin_no);

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
    } while (pMask[index] == true);
}

//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
void print_device_info()
{
    const int kb = 1024;
    const int mb = kb * kb;
    size_t free, total;
    

    int32_t device_count, cuda_version, driver_version;
    checkCudaErrors(cudaGetDeviceCount(&device_count));
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