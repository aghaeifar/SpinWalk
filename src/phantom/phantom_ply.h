
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_ply.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 20.01.2025
 * Descrip  : 
 * -------------------------------------------------------------------------- */

#ifndef PLY_H
#define PLY_H

#include "phantom_base.h"

namespace phantom
{
class ply : public phantom_base
{
    public:
    ply();
    ply(float fov_um, size_t resolution, float dChi, float Y, std::string filename_in, std::string filename_out = "shape.h5");
    virtual ~ply();

    virtual bool run(bool write_to_disk) override;
    virtual bool generate_mask_fieldmap();

    friend std::ostream& operator<<(std::ostream& os, const ply& obj);

    private: 
    std::vector<std::vector<float>> m_sphere_points;
    std::vector<float> m_sphere_radii;
    std::string filename_ply;
};

}

#endif // PLY_H