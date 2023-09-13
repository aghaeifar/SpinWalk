/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: kernels.cuh
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */


#ifndef _KERNELS_H_
#define _KERNELS_H_

#include <numeric>
#include "miscellaneous.h"

#define GAMMA  267515315. // rad/s.T

__global__ void cu_sim(const simulation_parameters *param, const float *pFieldMap, const bool *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint32_t spin_no = 0);

// scale position to mimic the different volume size
__global__ void cu_scalePos(float *scaled_xyz, float *initial_xyz, float scale, uint32_t size);

// generate random initial position
__global__ void cu_randPosGen(float *spin_position_xyz, simulation_parameters *param, bool *pMask, uint32_t spin_no = 0);


__host__  __device__ __forceinline__ uint64_t sub2ind(uint32_t x, uint32_t y, uint32_t z, uint32_t lenx, uint32_t leny)
{
    return (uint64_t((z-1)*lenx*leny) + (y-1)*lenx + x-1); // the last -1 is because of the C++ indexing starts from 0
}

void print_device_info();


template <typename T>
uint32_t is_masked(std::vector<T> &XYZ0, std::vector<char> &mask, simulation_parameters *param)
{
    float scale2grid[3];
    for(uint8_t i=0; i<3; i++)
        scale2grid[i] = (param->fieldmap_size[i]-1.) / param->sample_length[i];

    std::vector<uint8_t> mask_counter(XYZ0.size()/3, 0);
    #pragma omp parallel for
    for (int32_t i=0; i<XYZ0.size()/3; i++)
    {
        uint64_t index = sub2ind(ROUND(XYZ0[3*i] * scale2grid[0]), ROUND(XYZ0[3*i+1] * scale2grid[1]), ROUND(XYZ0[3*i+2] * scale2grid[2]), param->fieldmap_size[0], param->fieldmap_size[1]);
        if (index >= mask.size())
            std::cout << ERR_MSG << "Index out of range: " << index << " >= " << mask.size() << std::endl;
        mask_counter[i] = mask[index];
    }
    return std::accumulate(mask_counter.begin(), mask_counter.end(), 0);
}


#endif // _KERNELS_H_