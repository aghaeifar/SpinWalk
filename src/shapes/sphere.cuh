
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
    sphere(float fov_um, size_t resolution, float dChi, float Y, float radius_um = 50, size_t num_spheres = 5, std::string filename = "shape.h5");
    ~sphere();

    virtual bool run();
    virtual void set_sphere_parameters(float radius_um = 50, size_t num_spheres = 5);    
    virtual void generate_shapes();
    virtual void generate_mask_fieldmap();
    virtual void print_info();

    protected:

    private: 
        float *m_pSphere_points;
        float *m_pSphere_radius;
        float m_radius;
        size_t m_num_spheres;
};

#endif // SPHERE_H