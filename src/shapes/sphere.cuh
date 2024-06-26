
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: file_utils.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#ifndef SPHERE_H
#define SPHERE_H

#include "shape_base.cuh"

class sphere : public shape
{
    public:
    sphere();
    sphere(float fov_um, size_t resolution, float dChi, float Y, float radius_um = 50, float BVF = 10.0, bool is_seed_fixed=false, std::string filename = "shape.h5");
    ~sphere();

    virtual bool run();
    virtual void set_sphere_parameters(float radius_um = 50);    
    virtual void generate_shapes();
    virtual void generate_mask_fieldmap();
    virtual void print_info();

    protected: 

    private: 
    std::vector<float> m_sphere_points;
    std::vector<float> m_sphere_radii;
    float m_radius;
};

#endif // SPHERE_H