
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_base.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
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
        phantom_base(float fov_um, size_t resolution, float dChi, float Y, float BVF, int32_t seed, std::string filename);
        virtual ~phantom_base();
        void set_space(float fov_um, size_t resolution);
        void set_blood_parameters(float dChi, float Y, float BVF = 10.0); 
        void set_filename(std::string filename = "shape.h5");      
        virtual bool run(){return true;}; 
        virtual bool save() const;
        virtual bool create_grid();
        virtual void generate_shapes() = 0;
        virtual void generate_mask_fieldmap() = 0;

        friend std::ostream& operator<<(std::ostream& os, const phantom_base& obj);

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

   
float dot_product(const float *a, const float *b);
float norm(const float *a);
float norm_p2(const float *a);
void cross_product(const float *a, const float *b, float *c);
void normalize(float *a, float n = -1.0f);
void subtract(const float *a, const float *b, float *c);
void add(const float *a, const float *b, float *c);
void multiply(const float a, const float *b, float *c);
void multiply(const float *a, const float *b, float *c);
void copy(const float *a, float *b);
void roty(float theta, const float *m0, float *m1);

}

#endif // PHANTOM_BASE_H