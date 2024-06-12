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
#include "simulation_parameters.h"


#define GAMMA  267515315. // rad/s.T

__global__ void cu_sim(const simulation_parameters *param, const float *pFieldMap, const uint8_t *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint8_t *T);

// scale position to mimic the different volume size
__global__ void cu_scalePos(float *scaled_xyz, float *initial_xyz, float scale, uint64_t size);

// CUDA kernel to perform array multiplication with a constant
__global__ void cu_scaleArray(float *array, double scale, uint64_t size);

// generate random initial position
__global__ void cu_randPosGen(float *spin_position_xyz, simulation_parameters *param, const uint8_t *pMask, uint32_t spin_no = 0);

__host__  __device__ __forceinline__ int64_t sub2ind(int32_t x, int32_t y, int32_t z, int32_t lenx, int32_t leny)
{
    return (int64_t(z)*int64_t(lenx*leny) + y*lenx + x); 
}


#endif // _KERNELS_H_