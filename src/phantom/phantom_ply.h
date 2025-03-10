
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
#include <array>

namespace phantom
{
struct Vec3 {double x, y, z;
    Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
    Vec3 operator+(const double& other) const { return {x + other, y + other, z + other}; }
    Vec3 operator/(const double& other) const { return {x / other, y / other, z / other}; }
    Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
};
struct Triangle { Vec3 v0, v1, v2; };
// Axis-Aligned Bounding Box (AABB)
struct AABB { Vec3 min, max; 
    void expand(const Vec3& p) {
        min.x = std::min(min.x, p.x);
        min.y = std::min(min.y, p.y);
        min.z = std::min(min.z, p.z);
        max.x = std::max(max.x, p.x);
        max.y = std::max(max.y, p.y);
        max.z = std::max(max.z, p.z);
    }

    // Ray-AABB intersection. Valid only for ray_direction = {1., 0., 0.} do not change direction! 
    bool intersectsRay(const Vec3& origin) const {
        if(origin.y < min.y || origin.y > max.y || origin.z < min.z || origin.z > max.z)
            return false;
        return true;
    }
};

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
    AABB boundry;
    std::string filename_ply;
};

}

#endif // PLY_H