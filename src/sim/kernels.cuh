/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: kernels.cuh
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */


#ifndef KERNELS_H
#define KERNELS_H

#include <numeric>
#include "simulation_parameters.cuh"

#define ROW_MAJOR
// #define COL_MAJOR

namespace sim {
#ifdef __CUDACC__
__global__ 
#endif
void cu_sim(const parameters param, const parameters_uvec param_uvec, const float *pFieldMap, const uint8_t *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint8_t *T);

//
#ifdef __CUDACC__
__host__  __device__ 
#endif 
void sim(const parameters &param, const parameters_uvec &param_uvec, const float *pFieldMap, const uint8_t *pMask, const float *M0, const float *XYZ0, float *M1, float *XYZ1, uint8_t *T, uint32_t spin_no);


// data is stored row-major in the h5 file -> (x1,y1,z1); (x1,y1,z2); (x1,y1,z3)...(x1,y2,z1); (x1,y2,z2); (x1,y2,z3)...
#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
int64_t sub2ind(int64_t x, int64_t y, int64_t z, int64_t o, int64_t len_dim_x, int64_t len_dim_y, int64_t len_dim_z, int64_t len_dim_o)
{
#ifdef ROW_MAJOR
    return (x*len_dim_o*len_dim_z*len_dim_y + y*len_dim_z*len_dim_o + z*len_dim_o + o); 
#else
    return (o*len_dim_x*len_dim_y*len_dim_z + z*len_dim_x*len_dim_y + y*len_dim_x + x); // column-major
#endif
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
int64_t sub2ind(int64_t x, int64_t y, int64_t z, int64_t len_dim_x, int64_t len_dim_y, int64_t len_dim_z)
{
#ifdef ROW_MAJOR
    return (x*len_dim_z*len_dim_y + y*len_dim_z + z); 
#else
    return (z*len_dim_x*len_dim_y + y*len_dim_x + x); // column-major
#endif
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
int64_t sub2ind(int64_t x, int64_t y,int64_t len_dim_x, int64_t len_dim_y)
{
#ifdef ROW_MAJOR
    return (x*len_dim_y + y); 
#else
    return (y*len_dim_x + x); // column-major
#endif
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void xrot(float sin_theta, float cos_theta, const float *m0, float *m1)
{
    m1[0] = m0[0]; 
    m1[1] = cos_theta*m0[1] - sin_theta*m0[2];
    m1[2] = sin_theta*m0[1] + cos_theta*m0[2];
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void xrot(float theta, const float *m0, float *m1)
{
    float s = sinf(theta * DEG2RAD);
    float c = cosf(theta * DEG2RAD);
    xrot(s, c, m0, m1);
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void yrot(float sin_theta, float cos_theta, const float *m0, float *m1)
{
    m1[0] =  cos_theta*m0[0] + sin_theta*m0[2];
    m1[1] =  m0[1];
    m1[2] = -sin_theta*m0[0] + cos_theta*m0[2];
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void yrot(float theta, const float *m0, float *m1)
{
    float s = sinf(theta * DEG2RAD);
    float c = cosf(theta * DEG2RAD);
    yrot(s, c, m0, m1);
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void zrot(float sin_theta, float cos_theta, const float *m0, float *m1)
{
    m1[0] = cos_theta*m0[0] - sin_theta*m0[1];
    m1[1] = sin_theta*m0[0] + cos_theta*m0[1];
    m1[2] = m0[2];
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void zrot(float theta, const float *m0, float *m1)
{
    float s = sinf(theta * DEG2RAD);
    float c = cosf(theta * DEG2RAD);
    zrot(s, c, m0, m1);
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void xrot_withphase(float sin_theta, float cos_theta, float rf_phase, const float *m0, float *m1)
{
    if (rf_phase == 0.0f)
    {
        xrot(sin_theta, cos_theta, m0, m1);
        return;
    }
    if (rf_phase == 180.0)
    {
        xrot(-sin_theta, cos_theta, m0, m1);
        return;
    }
    if (rf_phase == 90.0)
    {
        yrot(sin_theta, cos_theta, m0, m1);
        return;
    }
    if (rf_phase == -90.0 || rf_phase == 270.0)
    {
        yrot(-sin_theta, cos_theta, m0, m1);
        return;
    }

    float m1_t[3];
    float s = sinf(rf_phase * DEG2RAD);
    float c = cosf(rf_phase * DEG2RAD);    
    zrot(-s, c, m0, m1);
    xrot(sin_theta, cos_theta, m1, m1_t);
    zrot(s, c, m1_t, m1);
    // MATLAB code:
    // fa = pi/5;
    // v = [1, 3, 0];
    // disp(axang2rotm([v, fa]))
    // P = atan2d(v(2), v(1));
    // disp(zrot(P) * xrot(rad2deg(fa)) * zrot(-P));
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void xrot_withphase(float theta, float rf_phase, const float *m0, float *m1)
{
    xrot_withphase(sin(theta*DEG2RAD), cos(theta*DEG2RAD), rf_phase, m0, m1);
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void relax(float e1, float e2, const float *m0, float *m1)
{
    m1[0] = m0[0] * e2;
    m1[1] = m0[1] * e2;
    m1[2] = 1. + e1*(m0[2] - 1.);
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void relax(float e1, float e2, float *m)
{
    relax(e1, e2, m, m);
}


} // namespace sim

#endif // _KERNELS_H_