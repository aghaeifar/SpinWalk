
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_ply.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 20.01.2022
 * Descrip  : 
 * -------------------------------------------------------------------------- */
#include <filesystem>
#include <random>
// boost includes
#include <boost/log/trivial.hpp> 

#include <highfive/highfive.hpp>
#include "barkeep.h"
#include "happly.h"
#include "phantom_ply.h"

// -------------------------------------------------------------------------- //
namespace phantom
{


bool is_inside_particle(const std::array<double, 3> &ray_origin, const std::vector<std::array<std::array<double, 3>, 3>> & triangles) {
    const double EPSILON = 1e-9;
    std::array<double ,3> ray_direction = {1., 0., 0.};

    // Möller–Trumbore algorithm for ray-triangle intersection: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    size_t num_intersection = 0;
    for(const auto &triangle : triangles){
        // Vertices of the triangle
        std::array<double, 3> v0 = triangle[0];
        std::array<double, 3> v1 = triangle[1];
        std::array<double, 3> v2 = triangle[2];

        // Find vectors for two edges sharing v0
        std::array<double, 3> edge1 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        std::array<double, 3> edge2 = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};

        // Compute determinant
        double ray_cross_e2[3];
        cross_product(ray_direction.data(), edge2.data(), ray_cross_e2);
        double det = dot_product(edge1.data(), ray_cross_e2);

        // If the determinant is near zero, the ray is parallel to the triangle
        if (std::fabs(det) < EPSILON) continue;

        // Calculate inverse determinant
        double inv_det= 1.0 / det;

        // Calculate the vector from v0 to the ray origin
        std::array<double, 3> s = {ray_origin[0] - v0[0], ray_origin[1] - v0[1], ray_origin[2] - v0[2]};

        // Calculate u parameter and test bounds
        double u = inv_det * dot_product(s.data(), ray_cross_e2);
        if (u < 0.0 || u > 1.0) continue;

        // Calculate v parameter and test bounds
        double s_cross_e1[3];
        cross_product(s.data(), edge1.data(), s_cross_e1);
        double v = inv_det * dot_product(ray_direction.data(), s_cross_e1);
        if (v < 0.0f || u + v > 1.0f) continue;

        // Calculate t to find out where the intersection point is on the line
        float t = inv_det * dot_product(edge2.data(), s_cross_e1);

        if (t > EPSILON)  // Ray intersection
            num_intersection++;
    }
    return num_intersection % 2;
}

ply::ply()
{

}

ply::ply(float fov_um, size_t resolution, float dChi, float Y, std::string filename_in, std::string filename_out)
: phantom_base(fov_um, resolution, dChi, Y, 0, 0, filename_out)
{
    filename_ply = filename_in;
}

ply::~ply()
{
}

bool ply::generate_mask_fieldmap()
{
    BOOST_LOG_TRIVIAL(info) << "Reading mesh..." << std::endl;
    
    happly::PLYData plyIn(filename_ply);

    // Get mesh-style data from the object
    std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
    std::vector<std::vector<size_t>> fInd   = plyIn.getFaceIndices<size_t>();
    for(const auto &v : fInd)
        if(v.size() != 3){
            BOOST_LOG_TRIVIAL(error) << "Only triangular mesh is supported!";
            return false;
        }
    
    std::vector<std::array<std::array<double, 3>, 3>> triangles;
    for(const auto &f : fInd)
        triangles.push_back({vPos[f[0]], vPos[f[1]], vPos[f[2]]});

    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;
    m_calc_fieldmap = false;
    m_fieldmap.resize(res3, 0);
    m_mask.resize(res3, 0);
    m_volume_fraction = 0;
    m_Y = -1;

    size_t c = 0;
    auto bar = barkeep::ProgressBar(&c, {.total = res1*res1, .message = "Progress", .style = barkeep::ProgressBarStyle::Rich,});
    
    for(int32_t pz=0; pz<res1; pz++)      
    for(int32_t py=0; py<res1; py++){
        c++;
    #pragma omp parallel for
    for(int32_t px=0; px<res1; px++)
    {
        std::array<double, 3> point;
        size_t p = px*res2 + py*res1 + pz;        
        point[0] = 1000*(m_grid[3*p] - m_fov/2.0); 
        point[1] = 1000*(m_grid[3*p+1] - m_fov/2.0); 
        point[2] = 1000*(m_grid[3*p+2] - m_fov/2.0);
        m_mask[p] = is_inside_particle(point, triangles);
    }
    }
    bar->done(); 

    return true;
}


std::ostream& operator<<(std::ostream& os, const ply& obj)
{
    os << static_cast<const phantom_base&>(obj);
    return os;
}

// -------------------------------------------------------------------------- //

bool ply::run(bool write_to_disk)
{
    BOOST_LOG_TRIVIAL(info) << *this << std::endl;
    if(create_grid() == false)
        return false;
    if(generate_mask_fieldmap() == false)
        return false;
    if(write_to_disk)
        if(save() == false)
            return false;
    return true;
}

}