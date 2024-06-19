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

#define ROW_MAJOR
// #define COL_MAJOR

#define GAMMA  267515315. // rad/s.T

__global__ void cu_sim(const simulation_parameters *param, const float *pFieldMap, const uint8_t *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint8_t *T);

// scale position to mimic the different volume size
__global__ void cu_scalePos(float *scaled_xyz, float *initial_xyz, float scale, uint64_t size);

// CUDA kernel to perform array multiplication with a constant
__global__ void cu_scaleArray(float *array, double scale, uint64_t size);

// generate random initial position
__global__ void cu_randPosGen(float *spin_position_xyz, simulation_parameters *param, const uint8_t *pMask, uint32_t spin_no = 0);

// data is stored row-major in the h5 file -> (x1,y1,z1); (x1,y1,z2); (x1,y1,z3)...(x1,y2,z1); (x1,y2,z2); (x1,y2,z3)...
__host__  __device__ __forceinline__ int64_t sub2ind(int64_t x, int64_t y, int64_t z, int64_t o, int64_t len_dim_x, int64_t len_dim_y, int64_t len_dim_z, int64_t len_dim_o)
{
#ifdef ROW_MAJOR
    return (x*len_dim_o*len_dim_z*len_dim_y + y*len_dim_z*len_dim_o + z*len_dim_o + o); 
#else
    return (o*len_dim_x*len_dim_y*len_dim_z + z*len_dim_x*len_dim_y + y*len_dim_x + x); // column-major
#endif
}

__host__  __device__ __forceinline__ int64_t sub2ind(int64_t x, int64_t y, int64_t z, int64_t len_dim_x, int64_t len_dim_y, int64_t len_dim_z)
{
#ifdef ROW_MAJOR
    return (x*len_dim_z*len_dim_y + y*len_dim_z + z); 
#else
    return (z*len_dim_x*len_dim_y + y*len_dim_x + x); // column-major
#endif
}

__host__  __device__ __forceinline__ int64_t sub2ind(int64_t x, int64_t y,int64_t len_dim_x, int64_t len_dim_y)
{
#ifdef ROW_MAJOR
    return (x*len_dim_y + y); 
#else
    return (y*len_dim_x + x); // column-major
#endif
}


#endif // _KERNELS_H_