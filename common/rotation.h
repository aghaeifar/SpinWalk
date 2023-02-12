
/* --------------------------------------------------------------------------
 * Project: 
 * File: 
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : 
 * -------------------------------------------------------------------------- */

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
void zrot(float theta, float *m0, float *m1)
{
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
void relax(float e1, float e2, float *m)
{
    m[0] *= e2;
    m[1] *= e2;
    m[2] = 1. + e1*(m[2] - 1.);
}