
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: file_utils.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#ifndef CYLINDER_H
#define CYLINDER_H

#include "shape_base.cuh"

class cylinder : public shape
{
    public:
    cylinder();
    cylinder(float fov_um, size_t resolution, float dChi, float Y, float radius_um = 50, size_t num_cylinders = 5, float orientation = -1.0f, std::string filename = "shape.h5");
    ~cylinder();

    virtual bool run();
    virtual void set_cylinder_parameters(float radius_um = 50, size_t num_cylinders = 5, float orientation = -1.0f);    
    virtual void generate_shapes();
    virtual void generate_mask_fieldmap();
    virtual void print_info();

    protected:

    private: 
        float *m_pCylinder_points;
        float *m_pCylinder_directions; // normalized vector
        float *m_pCylinder_radius;
        float m_radius, m_orientation;
        size_t m_num_cylinders;
};

#endif // CYLINDER_H