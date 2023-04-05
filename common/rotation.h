
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
    float s = sinf(theta);
    float c = cosf(theta);
    m1[0] = m0[0]; 
    m1[1] = c*m0[1] - s*m0[2];
    m1[2] = s*m0[1] + c*m0[2];
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__
#else
inline
#endif
void yrot(float theta, float *m0, float *m1)
{
    float s = sinf(theta);
    float c = cosf(theta);
    m1[0] = c*m0[0] + s*m0[2];
    m1[1] = m0[1];
    m1[2] = -s*m0[0] + c*m0[2];
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
    if(theta == 0.0f)
    {
        m1[0] = m0[0];
        m1[1] = m0[1];
        m1[2] = m0[2];
        return;
    }
    float s = sinf(theta);
    float c = cosf(theta);
    m1[0] = c*m0[0] - s*m0[1];
    m1[1] = s*m0[0] + c*m0[1];
    m1[2] = m0[2];
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__
#else
inline
#endif
void xrot_phasecycled(float sin_theta, float cos_theta, float phase_cycling, float *m0, float *m1)
{
    if (phase_cycling == 0.0f)
    {
        xrot(sin_theta, cos_theta, m0, m1);
        return;
    }

    if (phase_cycling == M_PI)
    {
        xrot(-sin_theta, cos_theta, m0, m1);
        return;
    }

    float m1_t[3];
    float s = sinf(phase_cycling);
    float c = cosf(phase_cycling);    
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