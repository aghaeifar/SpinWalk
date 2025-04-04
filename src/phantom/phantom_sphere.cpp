
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_sphere.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 15.05.2024
 * Descrip  : 
 * -------------------------------------------------------------------------- */
#include <filesystem>
#include <random>
// boost includes
#include <boost/log/trivial.hpp> 

#include <highfive/highfive.hpp>
#include "barkeep.h"
#include "phantom_sphere.h"

// -------------------------------------------------------------------------- //
namespace phantom
{

bool check_sphere_overlap(const std::vector<std::vector<float>>& m_sphere_points, const std::vector<float>& m_sphere_radii, const float* sph_pnt, float& radius, bool is_random_radius) {
    int c = 0;
    bool stop = false;

    #pragma omp parallel for shared(stop, radius)
    for (c = 0; c < m_sphere_radii.size(); c++) {
        if (stop) continue; // Skip computation if condition already met.

        float p2p1[3];
        subtract(sph_pnt, m_sphere_points[c].data(), p2p1);
        float distance = norm(p2p1);
        // Check the conditions
        if (distance <= m_sphere_radii[c] || distance <= radius) {
            #pragma omp critical
            {
                stop = true; // Signal all threads to stop further checks.
            }
        } else if (distance < m_sphere_radii[c] + radius) {
            if (!is_random_radius) {
                #pragma omp critical
                {
                    stop = true; // Signal all threads to stop further checks.
                }
            } else {
                #pragma omp critical
                {
                    if (!stop) 
                        radius = distance - m_sphere_radii[c]; // Adjust radius.
                }
            }
        }
    }
    return stop; // Return whether an overlap condition was found.
}

sphere::sphere()
{
    m_radius = 0;
}

sphere::sphere(float fov_um, size_t resolution, float dChi, float Y, float radius_um, float volume_fraction, int32_t seed, std::string filename)
: phantom_base(fov_um, resolution, dChi, Y, volume_fraction, seed, filename)
{
    set_sphere_parameters(radius_um);
}

sphere::~sphere()
{
}

void sphere::set_sphere_parameters(float radius)
{
    m_radius = radius;
}


bool sphere::generate_shapes()
{
     if(2*m_radius>=m_fov)
    {
        BOOST_LOG_TRIVIAL(error) << "Error: The radius of the cylinder is too large for the given FOV!\n";
        return false;
    }
    BOOST_LOG_TRIVIAL(info) << "Generating coordinates...for target BVF = " << m_volume_fraction << "% ...\n";    
    bool is_random_radius = m_radius < 0;
    float max_radius    = m_radius>0 ? m_radius:-m_radius;
    m_sphere_points.clear();
    m_sphere_radii.clear();
    float sph_pnt[3], radius ;

    std::minstd_rand gen(m_seed);
    std::uniform_real_distribution<float> dist(0.f, 1.f); 
      
    float distance, vol_sph = 0, vol_tol = m_fov*m_fov*m_fov;
    int32_t progress = 0;
    auto bar = barkeep::ProgressBar(&progress, {.total = 100, .message = "Generating (1/2)", .style = barkeep::ProgressBarStyle::Rich,});
    auto start = std::chrono::high_resolution_clock::now();
    while(progress < 100)
    {
        radius = is_random_radius ? dist(gen) * max_radius : max_radius;
        for (size_t i = 0; i < 3; i++) // generate a random point for a sphere which fit in the FOV
            sph_pnt[i] = dist(gen) * m_fov;        
        // check if sphere coordinate is ok
        if (check_sphere_overlap(m_sphere_points, m_sphere_radii, sph_pnt, radius, is_random_radius)) 
            continue;   

        vol_sph += 4*M_PI/3 * radius*radius*radius;
        m_sphere_points.push_back({sph_pnt[0], sph_pnt[1], sph_pnt[2]});
        m_sphere_radii.push_back(radius);  
        progress = 0.95 * 100*(100.*vol_sph/vol_tol/m_volume_fraction); // 0.95 is a factor to compensate for spheres in the boundary   
    }
    bar->done();

    auto end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) <<m_sphere_radii.size() << " coordinates generated successfully! Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << std::endl;
    return true;
}

bool sphere::generate_mask_fieldmap()
{
    BOOST_LOG_TRIVIAL(info) << "Generating spheres..." << std::endl;
    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;
    int32_t x_min, x_max, y_min, y_max, z_min, z_max;
    int32_t sph_rad_vox, sph_center_vox[3];

    BOOST_LOG_TRIVIAL(info) << "Allocating memory...";
    m_fieldmap.resize(m_calc_fieldmap ? res3:0, 0.f);
    m_mask.resize(res3, 0);
    BOOST_LOG_TRIVIAL(info) << "Done!\n";
    float v_size = m_fov / m_resolution;

    size_t c = 0;
    auto bar = barkeep::ProgressBar(&c, {.total = m_sphere_radii.size(), .message = "Generating (2/2)", .style = barkeep::ProgressBarStyle::Rich,});
    auto start = std::chrono::high_resolution_clock::now();
    for (c = 0; c < m_sphere_radii.size(); c++)
    {
        float *sph_center = m_sphere_points[c].data();
        float sph_rad   = m_sphere_radii[c];
        float sph_rad2  = sph_rad*sph_rad;
        sph_rad_vox = std::ceil(sph_rad / v_size)+1;
        sph_center_vox[0] = int32_t(sph_center[0]/v_size); 
        sph_center_vox[1] = int32_t(sph_center[1]/v_size);
        sph_center_vox[2] = int32_t(sph_center[2]/v_size);
        
        // find the bounding box of the sphere
        if (m_calc_fieldmap)
        {
            x_min = std::max(0, sph_center_vox[0] - sph_rad_vox*20);
            x_max = std::min((int32_t)m_resolution, sph_center_vox[0] + sph_rad_vox*20);
            y_min = std::max(0, sph_center_vox[1] - sph_rad_vox*20);
            y_max = std::min((int32_t)m_resolution, sph_center_vox[1] + sph_rad_vox*20);
            z_min = std::max(0, sph_center_vox[2] - sph_rad_vox*20);
            z_max = std::min((int32_t)m_resolution, sph_center_vox[2] + sph_rad_vox*20);

        } else 
        {
            x_min = std::max(0, sph_center_vox[0] - sph_rad_vox);
            x_max = std::min((int32_t)m_resolution, sph_center_vox[0] + sph_rad_vox + 2);
            y_min = std::max(0, sph_center_vox[1] - sph_rad_vox);
            y_max = std::min((int32_t)m_resolution, sph_center_vox[1] + sph_rad_vox + 2);
            z_min = std::max(0, sph_center_vox[2] - sph_rad_vox);
            z_max = std::min((int32_t)m_resolution, sph_center_vox[2] + sph_rad_vox + 2);
        }

        #pragma omp parallel for
        for(int32_t pz=z_min; pz<z_max; pz++)
        for(int32_t py=y_min; py<y_max; py++)
        for(int32_t px=x_min; px<x_max; px++)
        {
            size_t p = px*res2 + py*res1 + pz;
            float *grid = (float *)m_grid.data() + 3*p;
            float  p2p1[3], distance2, phi_c2, dp ;
            // distance between the points and vessel axis and vector from the projection point to the point
            subtract(grid, sph_center, p2p1);  // vector from the spatial points to the sphere center
            distance2 = norm_p2(p2p1);   // distance^2 between the points and vessel axis
            if (distance2 <=sph_rad2)
                m_mask[p] = 1;
            if (m_calc_fieldmap)
            {   // calculate the fieldmap from the sphere 
                // reference: https://doi.org/10.1002/nbm.1079
                dp = dot_product(p2p1, B0);   
                phi_c2 = dp * dp / distance2; // cos^2(phi)            
                m_fieldmap[p] += distance2>sph_rad2 ? 4*M_PI*(1-m_Y)*m_dChi * sph_rad2*sph_rad/distance2/sqrtf(distance2) * (phi_c2 - 1./3.) : 0.f;
            }
        }  
    } 
    bar->done();

    m_volume_fraction = std::accumulate(m_mask.begin(), m_mask.end(), 0) * 100.0 / m_mask.size();
    BOOST_LOG_TRIVIAL(info) << "Actual Volume Fraction = " << m_volume_fraction << "% ...\n";   
    auto end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Spheres generated successfully! " << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";
    return true;
}


std::ostream& operator<<(std::ostream& os, const sphere& obj)
{
    os << static_cast<const phantom_base&>(obj) 
       << "  Radius: " << obj.m_radius << " \xC2\xB5m\n"
       << "  Volume Fraction: " << obj.m_volume_fraction << "\n";
    return os;
}

// -------------------------------------------------------------------------- //

bool sphere::run(bool write_to_disk)
{
    BOOST_LOG_TRIVIAL(info) << *this << std::endl;
    if(create_grid() == false)
        return false;
    if(generate_shapes() == false)
        return false;
    if(generate_mask_fieldmap() == false)
        return false;
    if(write_to_disk)
        if(save() == false)
            return false;
    return true;
}

}