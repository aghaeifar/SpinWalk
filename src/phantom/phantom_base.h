
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_base.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : 
 * -------------------------------------------------------------------------- */

#ifndef PHANTOM_BASE_H
#define PHANTOM_BASE_H

#include <vector>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace phantom
{

class phantom_base
{
    public:
        phantom_base();
        phantom_base(float fov_um, size_t resolution, float dChi, float Y, float volume_fraction, int32_t seed, std::string filename);
        virtual ~phantom_base();
        void set_space(float fov_um, size_t resolution);
        void set_parameters(float dChi, float Y, float volume_fraction = 10.0); 
        void set_filename(std::string filename = "shape.h5");      
        virtual bool run(bool write_to_disk) = 0; 
        virtual bool save() const;
        virtual bool save(float fov[3], std::vector<size_t> resolution) const;
        virtual bool create_grid();
        float get_actual_volume_fraction() const {return m_volume_fraction;}

        friend std::ostream& operator<<(std::ostream& os, const phantom_base& obj);

    protected:
        std::vector<float> m_grid;
        std::vector<float> m_fieldmap;
        std::vector<int8_t> m_mask;
        size_t m_resolution;
        float m_fov, m_dChi, m_Y;
        float m_volume_fraction; // volume fraction
        std::string m_filename;
        float B0[3] = {0.f, 0.f, 1.f};
        bool m_calc_fieldmap;
        size_t m_seed;

    private: 
};


template< class T>
T dot_product(const T *a, const T *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

template< class T>
T norm(const T *a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

template< class T>
T norm_p2(const T *a)
{
    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]; // avoid sqrt
}

template< class T>
void cross_product(const T *a, const T *b, T *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

template< class T>
void normalize(T *a, T n=-1.0)
{
    if (n < 0)
        n = norm(a);
    a[0] /= n;
    a[1] /= n;
    a[2] /= n;
}

template< class T>
void subtract(const T *a, const T *b, T *c)
{
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

template< class T>
void add(const T *a, const T *b, T *c)
{
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

template< class T>
void multiply(const T a, const T *b, T *c)
{
    c[0] = a * b[0];
    c[1] = a * b[1];
    c[2] = a * b[2];
}

template< class T>
void multiply(const T *a, const T *b, T *c)
{
    c[0] = a[0] * b[0];
    c[1] = a[1] * b[1];
    c[2] = a[2] * b[2];
}

template< class T>
void copy(const T *a, T *b)
{
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
}

template< class T>
void roty(T theta, const T *m0, T *m1)
{
    T deg2rad = 0.0174532925199433; // = M_PI/180 
    T s = sin(theta * deg2rad);
    T c = cos(theta * deg2rad);
    m1[0] =  c*m0[0] + s*m0[2];
    m1[1] =  m0[1];
    m1[2] = -s*m0[0] + c*m0[2];
}

}

#endif // PHANTOM_BASE_H