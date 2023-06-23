
/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: rotation.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : 
 * -------------------------------------------------------------------------- */

#ifndef __ROTATION_H__
#define __ROTATION_H__

#include <cmath>

#define DEG2RAD 0.0174532925199433 // = M_PI/180 
#define RAD2DEG 57.2957795130823


#ifdef __CUDACC__
__host__  __device__ __forceinline__
#else
inline
#endif
void xrot(float sin_theta, float cos_theta, float *m0, float *m1)
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
void xrot(float theta, float *m0, float *m1)
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
void yrot(float sin_theta, float cos_theta, float *m0, float *m1)
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
void yrot(float theta, float *m0, float *m1)
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
void zrot(float sin_theta, float cos_theta, float *m0, float *m1)
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
void zrot(float theta, float *m0, float *m1)
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
void xrot_withphase(float sin_theta, float cos_theta, float rf_phase, float *m0, float *m1)
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
void xrot_withphase(float theta, float rf_phase, float *m0, float *m1)
{
    xrot_withphase(sin(theta*DEG2RAD), cos(theta*DEG2RAD), rf_phase, m0, m1);
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__
#else
inline
#endif
void relax(float e1, float e2, float *m0, float *m1)
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

#endif // __ROTATION_H__