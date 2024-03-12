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

__global__ void cu_sim(const simulation_parameters *param, const float *pFieldMap, const uint8_t *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1);

// scale position to mimic the different volume size
__global__ void cu_scalePos(float *scaled_xyz, float *initial_xyz, float scale, uint64_t size);

// CUDA kernel to perform array multiplication with a constant
__global__ void cu_scaleArray(float *array, float scale, uint64_t size);

// generate random initial position
__global__ void cu_randPosGen(float *spin_position_xyz, simulation_parameters *param, const uint8_t *pMask, uint32_t spin_no = 0);

__host__  __device__ __forceinline__ uint64_t sub2ind(uint32_t x, uint32_t y, uint32_t z, uint32_t lenx, uint32_t leny)
{
    return (uint64_t((z-1)*lenx*leny) + (y-1)*lenx + x-1); // the last -1 is because of the C++ indexing starts from 0
}

uint32_t getDeviceCount();
void print_device_info();
bool check_memory_size(size_t required_size_MB);

template <typename T>
bool is_masked(std::vector<T> &XYZ0, std::vector<uint8_t> &mask, simulation_parameters *param)
{
    float scale2grid[3];
    for(uint8_t i=0; i<3; i++)
        scale2grid[i] = (param->fieldmap_size[i]-1.) / param->sample_length[i];

    for (int32_t i=0; i<XYZ0.size()/3; i++)
    {
        uint64_t index = sub2ind(ROUND(XYZ0[3*i] * scale2grid[0]), ROUND(XYZ0[3*i+1] * scale2grid[1]), ROUND(XYZ0[3*i+2] * scale2grid[2]), param->fieldmap_size[0], param->fieldmap_size[1]);
        if (index >= mask.size())
        {
            std::cout << ERR_MSG << "For "<<i<< "th element, index is out of range: " << index << " >= " << mask.size() << std::endl;
            return false;
        }
        if (mask[index] != 0 && param->enMultiTissue == false)
        {
            std::cout << ERR_MSG << " " << i<< "th element located in the masked!" << std::endl;
            return false;
        }
    }
    return true;
}


#endif // _KERNELS_H_