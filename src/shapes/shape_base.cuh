
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
        shape(float fov_um, size_t resolution, float dChi, float Y, float BVF, bool is_seed_fixed, std::string filename);
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
        bool m_random_seed;

    private: 
};

#endif // SHAPE_BASE_H