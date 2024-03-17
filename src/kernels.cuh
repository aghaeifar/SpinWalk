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


#endif // _KERNELS_H_