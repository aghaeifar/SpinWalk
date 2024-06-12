
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
    m_num_spheres = 0;
    m_pSphere_points = nullptr;
    m_pSphere_radius = nullptr;
}

sphere::sphere(float fov_um, size_t resolution, float dChi, float Y, float radius_um, size_t num_spheres, std::string filename)
: shape(fov_um, resolution, dChi, Y, filename)
{
    m_pSphere_points = nullptr;
    m_pSphere_radius = nullptr;
    set_sphere_parameters(radius_um, num_spheres);
}

sphere::~sphere()
{
    if (m_pSphere_points != nullptr)
        delete[] m_pSphere_points;
    if (m_pSphere_radius != nullptr)
        delete[] m_pSphere_radius;
}

void sphere::set_sphere_parameters(float radius, size_t num_spheres)
{
    m_radius = radius;
    m_num_spheres = num_spheres;
}

void sphere::generate_shapes()
{
    std::cout << "Generating spheres...\n";
    std::cout << std::fixed << std::setprecision(5);
    if (m_pSphere_points != nullptr)
        delete[] m_pSphere_points;
    if (m_pSphere_radius != nullptr)
        delete[] m_pSphere_radius; 

    bool is_random_radius = m_radius < 0;
    float max_radius    = m_radius>0 ? m_radius:-m_radius;
    m_pSphere_points      = new float[m_num_spheres * 3];
    m_pSphere_radius      = new float[m_num_spheres];
    float *cyl_pnt ;

    // srandom engine
    std::random_device rd; // Seed
    std::mt19937 gen(rd()); // Mersenne Twister generator
    std::uniform_real_distribution<> dist(0.f, 1.f); 
      
    size_t sph_counter = 0;
    tqdm bar;
    auto start = std::chrono::high_resolution_clock::now();
    while(sph_counter < m_num_spheres)
    {
        cyl_pnt = m_pSphere_points     + 3*sph_counter;
        // generate a random point and direction
        for (size_t i = 0; i < 3; i++)
            cyl_pnt[i] = dist(gen) * m_fov;
  
        m_pSphere_radius[sph_counter] = is_random_radius ? dist(gen) * max_radius : max_radius;

        // check if spheres overlap
        std::vector<float> cd(m_num_spheres, 0.f); 
        // #pragma omp parallel for // parallel loop is slower than serial loop! why?
        for (size_t c=0; c<sph_counter; c++)
        {   
            float p2p1[3];
            subtract(cyl_pnt, m_pSphere_points + 3*c, p2p1);
            cd[c] = norm(p2p1) - (m_pSphere_radius[sph_counter] + m_pSphere_radius[c]);
        }
        // all false
        if (std::none_of(cd.begin(), cd.end(), [](float v) { return v<0; }) == false)
            continue;
        
        bar.progress(sph_counter, m_num_spheres);
        sph_counter++;
    }
    bar.finish();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "spheres generated successfully! Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << std::endl;
}

void sphere::generate_mask_fieldmap()
{
    std::cout << "Generating fieldmap..." << std::endl;
    if (m_pFieldmap != nullptr)
        delete[] m_pFieldmap;
    if (m_pMask != nullptr)
        delete[] m_pMask;
    
    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;

    std::cout<<"Alocating memory..."<<std::endl;
    m_pFieldmap = new float[res3];
    m_pMask = new uint8_t[res3];
    std::fill(m_pFieldmap, m_pFieldmap + res3, 0.f);
    std::fill(m_pMask, m_pMask + res3, 0);

    tqdm bar;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t c = 0; c < m_num_spheres; c++)
    {
        float *sph_center = m_pSphere_points + 3*c;
        float sph_rad  = m_pSphere_radius[c];
        #pragma omp parallel for
        for(size_t p=0; p<res3; p++)
        {
            float *grid = m_pGrid + 3*p;
            float  p2p1[3], distance, phi_c ;
            // distance between the points and vessel axis and vector from the projection point to the point
            subtract(grid, sph_center, p2p1);  // vector from the spatial points to the sphere center
            distance = norm(p2p1);             // distance between the points and vessel axis
            phi_c = dot_product(p2p1, B0) / distance; // cos(phi)
            // calculate the fieldmap from the vessel 
            distance = sph_rad / distance; 
            m_pFieldmap[p] += distance<=1 ? 4*M_PI*(1-m_Y)*m_dChi*distance*distance*distance*(phi_c*phi_c - 1./3.) : 0.f;
            m_pMask[p] = distance<=1.f ? m_pMask[p] : 1;
        }        
        bar.progress(c, m_num_spheres);
    } 
    bar.finish();
    m_BVF = std::accumulate(m_pMask, m_pMask+res3, 0) * 100.0 / res3;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Fieldmaps generated successfully! BVF: " << m_BVF << "% Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";

}

void sphere::print_info()
{
    shape::print_info();
    std::cout << "  Radius: " << m_radius << " um\n";
    std::cout << "  Number of spheres: " << m_num_spheres << "\n";
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