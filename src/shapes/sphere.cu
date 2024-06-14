
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
#include "tqdm.h"
#include "sphere.cuh"
#include "basic_functions.cuh"

#define epsilon 1e-6

// -------------------------------------------------------------------------- //
sphere::sphere()
{
    m_radius = 0;
    m_BVF = 10.0;
}

sphere::sphere(float fov_um, size_t resolution, float dChi, float Y, float radius_um, float BVF, std::string filename)
: shape(fov_um, resolution, dChi, Y, filename)
{
    set_sphere_parameters(radius_um, BVF);
}

sphere::~sphere()
{
}

void sphere::set_sphere_parameters(float radius, float BVF)
{
    m_radius = radius;
    m_BVF = BVF;
}

void sphere::generate_shapes()
{
    std::cout << "Generating coordinates for BVF = " << m_BVF << "% ...\n";
    std::cout << std::fixed << std::setprecision(5);
    
    bool is_random_radius = m_radius < 0;
    float max_radius    = m_radius>0 ? m_radius:-m_radius;
    m_sphere_points.clear();
    m_sphere_radii.clear();
    float cyl_pnt[3], radius ;

    // srandom engine
    std::random_device rd; // Seed
    std::mt19937 gen(rd()); // Mersenne Twister generator
    std::uniform_real_distribution<> dist(0.f, 1.f); 
      
    float distance, vol_sph = 0, vol_tol = m_fov*m_fov*m_fov;
    int32_t progress = 0;
    tqdm bar;
    auto start = std::chrono::high_resolution_clock::now();
    while(progress < 100)
    {
        start_label:
        radius = is_random_radius ? dist(gen) * max_radius : max_radius;
        for (size_t i = 0; i < 3; i++) // generate a random point for a sphere which fit in the FOV
            cyl_pnt[i] = dist(gen) * (m_fov-2*radius) + radius;        
        // check if sphere coordinate is ok
        for (size_t c=0; c<m_sphere_radii.size(); c++)
        {   
            float p2p1[3];
            subtract(cyl_pnt, &m_sphere_points[3*c], p2p1);
            distance = norm(p2p1);
            // if the sphere is inside another sphere, generate a new sphere
            if(distance <= m_sphere_radii[c] ||  distance <= radius)
                goto start_label; // I hate to use goto statement, but let's do it for now because I want to combine break/continue for two loops and start from the beginning
            // adjust the radius of the sphere to avoid overlap
            if (distance < m_sphere_radii[c] + radius)
                radius = distance - m_sphere_radii[c];
        }

        vol_sph += 4*M_PI/3 * radius*radius*radius;
        m_sphere_points.insert(m_sphere_points.end(), cyl_pnt, cyl_pnt+3);
        m_sphere_radii.push_back(radius);  
        progress = 100*(100.*vol_sph/vol_tol/m_BVF);     
        bar.progress(progress, 100);
    }
    bar.finish();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout <<m_sphere_radii.size() << " coordinates generated successfully! Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << std::endl;

}

void sphere::generate_mask_fieldmap()
{
    std::cout << "Generating spheres..." << std::endl;

    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;

    std::cout<<"Allocating memory..."<<std::endl;
    m_fieldmap.resize(m_calc_fieldmap ? res3:0, 0.f);
    m_mask.resize(res3, 0);

    tqdm bar;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t c = 0; c < m_sphere_radii.size(); c++)
    {
        float *sph_center = &m_sphere_points[3*c];
        float sph_rad  = m_sphere_radii[c];
        #pragma omp parallel for
        for(size_t p=0; p<res3; p++)
        {
            float *grid = (float *)m_grid.data() + 3*p;
            float  p2p1[3], distance, phi_c ;
            // distance between the points and vessel axis and vector from the projection point to the point
            subtract(grid, sph_center, p2p1);  // vector from the spatial points to the sphere center
            distance = norm(p2p1);   // distance between the points and vessel axis
            if (distance <=sph_rad)
                m_mask[p] = 1;
            if (m_calc_fieldmap) 
            {   // calculate the fieldmap from the vessel 
                phi_c = dot_product(p2p1, B0) / distance; // cos(phi)            
                m_fieldmap[p] += distance<=1 ? 4*M_PI*(1-m_Y)*m_dChi * sph_rad*sph_rad*sph_rad/distance/distance/distance * (phi_c*phi_c - 1./3.) : 0.f;
            }
        }        
        bar.progress(c, m_sphere_radii.size());
    } 
    bar.finish();
    m_BVF = std::accumulate(m_mask.begin(), m_mask.end(), 0) * 100.0 / res3;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Spheres generated successfully! BVF: " << m_BVF << "% Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";

}

void sphere::print_info()
{
    shape::print_info();
    std::cout << "  Radius: " << m_radius << " um\n";
    std::cout << "  BVF " << m_BVF << "\n";
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