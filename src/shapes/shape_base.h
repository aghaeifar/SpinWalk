
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: file_utils.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#ifndef SHAPE_BASE_H
#define SHAPE_BASE_H

#include <vector>

class shape
{
    public:
        shape();
        shape(float fov_um, size_t resolution, float dChi, float Y, float BVF, int32_t seed, std::string filename);
        ~shape();
        void set_space(float fov_um, size_t resolution);
        void set_blood_parameters(float dChi, float Y, float BVF = 10.0); 
        void set_filename(std::string filename = "shape.h5");      
        virtual bool run(){return true;}; 
        virtual bool save();
        virtual bool create_grid();
        virtual void generate_shapes() = 0;
        virtual void generate_mask_fieldmap() = 0;
        virtual void print_info();

    protected:
        std::vector<float> m_grid;
        std::vector<float> m_fieldmap;
        std::vector<int8_t> m_mask;
        size_t m_resolution;
        float m_fov, m_dChi, m_Y;
        float m_BVF; // blood volume fraction
        std::string m_filename;
        float B0[3] = {0.f, 0.f, 1.f};
        bool m_calc_fieldmap;
        size_t m_seed;

    private: 
};


namespace shapes_functions
{      
    inline float dot_product(const float *a, const float *b)
    {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }

    inline float norm(const float *a)
    {
        return sqrtf(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    }

    inline float norm_p2(const float *a)
    {
        return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]; // avoid sqrt
    }

    inline void cross_product(const float *a, const float *b, float *c)
    {
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
    }

    inline void normalize(float *a, float n = -1.0f)
    {
        if (n < 0)
            n = norm(a);
        a[0] /= n;
        a[1] /= n;
        a[2] /= n;
    }

    inline void subtract(const float *a, const float *b, float *c)
    {
        c[0] = a[0] - b[0];
        c[1] = a[1] - b[1];
        c[2] = a[2] - b[2];
    }

    inline void add(const float *a, const float *b, float *c)
    {
        c[0] = a[0] + b[0];
        c[1] = a[1] + b[1];
        c[2] = a[2] + b[2];
    }

    inline void multiply(const float a, const float *b, float *c)
    {
        c[0] = a * b[0];
        c[1] = a * b[1];
        c[2] = a * b[2];
    }

    inline void multiply(const float *a, const float *b, float *c)
    {
        c[0] = a[0] * b[0];
        c[1] = a[1] * b[1];
        c[2] = a[2] * b[2];
    }

    inline void copy(const float *a, float *b)
    {
        b[0] = a[0];
        b[1] = a[1];
        b[2] = a[2];
    }

    inline void roty(float theta, const float *m0, float *m1)
    {
        float deg2rad = 0.0174532925199433; // = M_PI/180 
        float s = sinf(theta * deg2rad);
        float c = cosf(theta * deg2rad);
        m1[0] =  c*m0[0] + s*m0[2];
        m1[1] =  m0[1];
        m1[2] = -s*m0[0] + c*m0[2];
    }
}

#endif // SHAPE_BASE_H