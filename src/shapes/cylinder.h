
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

#include "shape_base.h"

class cylinder : public shape
{
    public:
    cylinder();
    cylinder(float fov_um, size_t resolution, float dChi, float Y, float radius_um = 50, float BVF = 10., float orientation = -1.0f, int32_t seed=-1, std::string filename = "shape.h5");
    ~cylinder();

    virtual bool run();
    virtual void set_cylinder_parameters(float radius_um = 50, float orientation = -1.0f);    
    virtual void generate_shapes();
    virtual void generate_mask_fieldmap();
    virtual void print_info();
    virtual float calculate_volume(float *cyl_pnt, float radius);

    protected:

    private: 
    std::vector<std::vector<float>> m_cylinder_points;
    std::vector<float> m_cylinder_radii;
    float m_radius, m_orientation;
};

#endif // CYLINDER_H