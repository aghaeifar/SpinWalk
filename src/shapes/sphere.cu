
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: arbitrary_gradient.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 15.05.2024
 * Descrip  : 
 * -------------------------------------------------------------------------- */
#include <highfive/highfive.hpp>
#include <filesystem>
#include <random>
#include "indicators.hpp"
#include "sphere.cuh"
#include "basic_functions.cuh"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace indicators;
// -------------------------------------------------------------------------- //
sphere::sphere()
{
    m_radius = 0;
}

sphere::sphere(float fov_um, size_t resolution, float dChi, float Y, float radius_um, float BVF, int32_t seed, std::string filename)
: shape(fov_um, resolution, dChi, Y, BVF, seed, filename)
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

void sphere::generate_shapes()
{
     if(2*m_radius>=m_fov)
    {
        std::cerr << "Error: The radius of the cylinder is too large for the given FOV!\n";
        return;
    }
    std::cout << "Generating coordinates...for target BVF = " << m_BVF << "% ...\n";    
    bool is_random_radius = m_radius < 0;
    float max_radius    = m_radius>0 ? m_radius:-m_radius;
    m_sphere_points.clear();
    m_sphere_radii.clear();
    float sph_pnt[3], radius ;

    // srandom engine
    std::mt19937 gen(m_seed); // Mersenne Twister generator
    std::uniform_real_distribution<float> dist(0.f, 1.f); 
      
    float distance, vol_sph = 0, vol_tol = m_fov*m_fov*m_fov;
    int32_t progress = 0;
    ProgressBar bar{option::ShowPercentage{true}, option::Start{"["}, option::Fill{"="}, option::Lead{">"}, option::End{"]"}};
    auto start = std::chrono::high_resolution_clock::now();
    while(progress < 100)
    {
        radius = is_random_radius ? dist(gen) * max_radius : max_radius;
        for (size_t i = 0; i < 3; i++) // generate a random point for a sphere which fit in the FOV
            sph_pnt[i] = dist(gen) * m_fov;        
        // check if sphere coordinate is ok
        size_t c;
        for (c=0; c<m_sphere_radii.size(); c++)
        {   
            float p2p1[3];
            subtract(sph_pnt, &m_sphere_points[3*c], p2p1);
            distance = norm(p2p1);
            // if the sphere is inside another sphere, generate a new sphere
            if(distance <= m_sphere_radii[c] ||  distance <= radius)
                break;
            // adjust the radius of the sphere to avoid overlap
            if (distance < m_sphere_radii[c] + radius)
            {
                if (!is_random_radius)
                    break;            
                radius = distance - m_sphere_radii[c];
            }
        }
        if (c < m_sphere_radii.size())
            continue;

        vol_sph += 4*M_PI/3 * radius*radius*radius;
        m_sphere_points.insert(m_sphere_points.end(), sph_pnt, sph_pnt+3);
        m_sphere_radii.push_back(radius);  
        progress = 0.95 * 100*(100.*vol_sph/vol_tol/m_BVF); // 0.95 is a factor to compensate for spheres in the boundary   
        bar.set_progress(progress);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout <<m_sphere_radii.size() << " coordinates generated successfully! Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << std::endl;

}

void sphere::generate_mask_fieldmap()
{
    std::cout << "Generating spheres..." << std::endl;
    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;
    int32_t x_min, x_max, y_min, y_max, z_min, z_max;
    int32_t sph_rad_vox, sph_center_vox[3];

    std::cout<<"Allocating memory...";
    m_fieldmap.resize(m_calc_fieldmap ? res3:0, 0.f);
    m_mask.resize(res3, 0);
    std::cout<<"Done!\n";
    float v_size = m_fov / m_resolution;

    std::cout << "Generating...\n";
    ProgressBar bar{option::ShowPercentage{true}, option::Start{"["}, option::Fill{"="}, option::Lead{">"}, option::End{"]"}, option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t c = 0; c < m_sphere_radii.size(); c++)
    {
        float *sph_center = &m_sphere_points[3*c];
        float sph_rad   = m_sphere_radii[c];
        float sph_rad2  = sph_rad*sph_rad;
        sph_rad_vox = std::ceil(sph_rad / v_size)+1;
        sph_center_vox[0] = int32_t(sph_center[0]/v_size); 
        sph_center_vox[1] = int32_t(sph_center[1]/v_size);
        sph_center_vox[2] = int32_t(sph_center[2]/v_size);
        
        // find the bounding box of the sphere
        if (m_calc_fieldmap)
        {
            x_min = std::max(0, sph_center_vox[0] - sph_rad_vox*4);
            x_max = std::min((int32_t)m_resolution, sph_center_vox[0] + sph_rad_vox*4);
            y_min = std::max(0, sph_center_vox[1] - sph_rad_vox*4);
            y_max = std::min((int32_t)m_resolution, sph_center_vox[1] + sph_rad_vox*4);
            z_min = std::max(0, sph_center_vox[2] - sph_rad_vox*4);
            z_max = std::min((int32_t)m_resolution, sph_center_vox[2] + sph_rad_vox*4);

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
            {   // calculate the fieldmap from the vessel 
                dp = dot_product(p2p1, B0);   
                phi_c2 = dp * dp / distance2; // cos^2(phi)            
                m_fieldmap[p] += distance2>sph_rad2 ? 4*M_PI*(1-m_Y)*m_dChi * sph_rad2*sph_rad/distance2/sqrtf(distance2) * (phi_c2 - 1./3.) : 0.f;
            }
        }  
        bar.set_progress(100 * (c+1)/float(m_sphere_radii.size()));
    } 

    m_BVF = std::accumulate(m_mask.begin(), m_mask.end(), 0) * 100.0 / m_mask.size();
    std::cout << "Actual Volume Fraction = " << m_BVF << "% ...\n";   
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Spheres generated successfully! " << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";
    
}

void sphere::print_info()
{
    shape::print_info();
    std::cout << "  Radius: " << m_radius << " um\n";
    std::cout << "  Volume Fraction " << m_BVF << "\n";
}

// -------------------------------------------------------------------------- //

bool sphere::run()
{
    print_info();
    create_grid();
    generate_shapes();
    generate_mask_fieldmap();
    save();
    return true;
}