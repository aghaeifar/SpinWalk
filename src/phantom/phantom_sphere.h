
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_sphere.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : 
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
    sphere(float fov_um, size_t resolution, float dChi, float Y, float radius_um = 50, float volume_fraction = 10.0, int32_t seed = -1, std::string filename = "shape.h5");
    virtual ~sphere();

    virtual bool run(bool write_to_disk) override;
    virtual void set_sphere_parameters(float radius_um = 50);    
    virtual bool generate_shapes();
    virtual bool generate_mask_fieldmap();

    friend std::ostream& operator<<(std::ostream& os, const sphere& obj);

    private: 
    std::vector<std::vector<float>> m_sphere_points;
    std::vector<float> m_sphere_radii;
    float m_radius;
};

}

#endif // SPHERE_H