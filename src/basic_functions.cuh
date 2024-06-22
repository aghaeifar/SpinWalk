
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: rotation.cuh
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


#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
float dot_product(const float *a, const float *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
float norm(const float *a)
{
    return sqrtf(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
float norm_p2(const float *a)
{
    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]; // avoid sqrt
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void cross_product(const float *a, const float *b, float *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void normalize(float *a, float n = -1.0f)
{
    if (n < 0)
        n = norm(a);
    a[0] /= n;
    a[1] /= n;
    a[2] /= n;
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void subtract(const float *a, const float *b, float *c)
{
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void add(const float *a, const float *b, float *c)
{
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void multiply(const float a, const float *b, float *c)
{
    c[0] = a * b[0];
    c[1] = a * b[1];
    c[2] = a * b[2];
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void multiply(const float *a, const float *b, float *c)
{
    c[0] = a[0] * b[0];
    c[1] = a[1] * b[1];
    c[2] = a[2] * b[2];
}

#ifdef __CUDACC__
__host__  __device__ __forceinline__ 
#else
inline
#endif
void copy(const float *a, float *b)
{
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
}

#endif // __ROTATION_H__