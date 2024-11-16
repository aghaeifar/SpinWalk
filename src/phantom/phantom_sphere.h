
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

#include "phantom_base.h"

namespace phantom
{
class sphere : public phantom_base
{
    public:
    sphere();
    sphere(float fov_um, size_t resolution, float dChi, float Y, float radius_um = 50, float BVF = 10.0, int32_t seed = -1, std::string filename = "shape.h5");
    virtual ~sphere();

    virtual bool run() override;
    virtual void set_sphere_parameters(float radius_um = 50);    
    virtual void generate_shapes();
    virtual void generate_mask_fieldmap();

    friend std::ostream& operator<<(std::ostream& os, const sphere& obj);

    private: 
    std::vector<float> m_sphere_points;
    std::vector<float> m_sphere_radii;
    float m_radius;
};

}

#endif // SPHERE_H