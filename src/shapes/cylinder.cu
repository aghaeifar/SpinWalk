
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
#include "cylinder.cuh"
#include "basic_functions.cuh"

#define epsilon 1e-6

cylinder::cylinder()
{
    m_radius = 0;
    m_orientation = 0;
    m_num_cylinders = 0;
    m_pCylinder_points = nullptr;
    m_pCylinder_directions = nullptr;
    m_pCylinder_radius = nullptr;
}

cylinder::cylinder(float fov_um, size_t resolution, float dChi, float Y, float radius_um, size_t num_cylinders, float orientation, std::string filename)
: shape(fov_um, resolution, dChi, Y, filename)
{
    m_pCylinder_points = nullptr;
    m_pCylinder_directions = nullptr;
    m_pCylinder_radius = nullptr;
    set_cylinder_parameters(radius_um, num_cylinders, orientation);
}

cylinder::~cylinder()
{
    if (m_pCylinder_points != nullptr)
        delete[] m_pCylinder_points;
    if (m_pCylinder_directions != nullptr)
        delete[] m_pCylinder_directions;
    if (m_pCylinder_radius != nullptr)
        delete[] m_pCylinder_radius;
}

void cylinder::set_cylinder_parameters(float radius, size_t num_cylinders, float orientation)
{
    m_radius = radius;
    m_num_cylinders = num_cylinders;
    m_orientation = orientation;
}

void cylinder::generate_shapes()
{
    std::cout << "Generating cylinders...\n";
    std::cout << std::fixed << std::setprecision(5);
    if (m_pCylinder_directions != nullptr)
        delete[] m_pCylinder_directions;
    if (m_pCylinder_points != nullptr)
        delete[] m_pCylinder_points;
    if (m_pCylinder_radius != nullptr)
        delete[] m_pCylinder_radius; 

    bool is_random_orientation = m_orientation < 0;
    bool is_random_radius = m_radius < 0;
    float max_radius    = m_radius>0 ? m_radius:-m_radius;
    m_pCylinder_directions  = new float[m_num_cylinders * 3];
    m_pCylinder_points      = new float[m_num_cylinders * 3];
    m_pCylinder_radius      = new float[m_num_cylinders];
    float *cyl_dir, *cyl_pnt ;

    // set the cylinders orientation if they are parallel
    float d[3];
    xrot(m_orientation, B0, d); // a vector with the specified angle with respect to vector B0

    // srandom engine
    std::random_device rd; // Seed
    std::mt19937 gen(rd()); // Mersenne Twister generator
    // std::mt19937 gen(0); // Mersenne Twister generator
    std::uniform_real_distribution<> dist(0.f, 1.f); 
      
    size_t cyl_counter = 0;
    tqdm bar;
    auto start = std::chrono::high_resolution_clock::now();
    while(cyl_counter < m_num_cylinders)
    {
        cyl_pnt = m_pCylinder_points     + 3*cyl_counter;
        cyl_dir = m_pCylinder_directions + 3*cyl_counter;
        // generate a random point and direction
        for (size_t i = 0; i < 3; i++)
        {
            cyl_pnt[i] = dist(gen) * m_fov;
            cyl_dir[i] = dist(gen) - 0.5;
        }
        m_pCylinder_radius[cyl_counter] = is_random_radius ? dist(gen) * max_radius : max_radius;
        normalize(cyl_dir);

        // check if cylinders overlap
        std::vector<float> cd(m_num_cylinders, 0.f); 
        // #pragma omp parallel for // parallel loop is slower than serial loop! why?
        for (size_t c=0; c<cyl_counter; c++)
        {   
            float p2p1[3], d2d1[3], d2d1_norm, distance;
            subtract(cyl_pnt, m_pCylinder_points + 3*c, p2p1);
            cross_product(cyl_dir, m_pCylinder_directions + 3*c, d2d1);
            d2d1_norm = norm(d2d1);
            if (d2d1_norm < epsilon) // cylinders are parallel
            {
                float temp[3];
                cross_product(p2p1, cyl_dir, temp);
                distance = norm(temp); // norm(temp) / curr_cyl_norm   where  curr_cyl_norm= 1.0;
            }
            else
            {
                distance = dot_product(p2p1, d2d1) / d2d1_norm;
                distance = distance<0 ? -distance : distance; // abs
            }
            cd[c] = distance - (m_pCylinder_radius[cyl_counter] + m_pCylinder_radius[c]);
        }
        // all false
        if (std::none_of(cd.begin(), cd.end(), [](float v) { return v<0; }) == false)
            continue;
        
        bar.progress(cyl_counter, m_num_cylinders);
        cyl_counter++;
    }
    bar.finish();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Cylinders generated successfully! Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << std::endl;
}

void cylinder::generate_mask_fieldmap()
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
    float theta_c, theta_s2;
    // project B0 from the projection point to the point to the plane perpendicular to the vessel axis
    float B0_prj[3];    
    tqdm bar;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t c = 0; c < m_num_cylinders; c++)
    {
        float *cyl_dir = m_pCylinder_directions + 3*c;
        float *cyl_pnt = m_pCylinder_points + 3*c;
        float cyl_rad  = m_pCylinder_radius[c];
        // project B0 from the projection point to the point to the plane perpendicular to the vessel axis  
        theta_c  = dot_product(cyl_dir, B0); // angle between axis of vessel and B0 (in radian)
        theta_s2 = 1.0 - theta_c * theta_c;
        multiply(theta_c, cyl_dir, B0_prj);
        subtract(B0, B0_prj, B0_prj);
        normalize(B0_prj);

        #pragma omp parallel for
        for(size_t p=0; p<res3; p++)
        {
            float *grid = m_pGrid + 3*p;
            float temp[3], perpendicular[3], distance, phi_c, phi_2c2_1 ;
            // distance between the points and vessel axis and vector from the projection point to the point
            subtract(grid, cyl_pnt, temp);  // vector from the spatial points to the cylinder point
            multiply(dot_product(cyl_dir, temp), cyl_dir, temp); // project vector temp onto the cylinder direction vector
            add(temp, cyl_pnt, temp);               // projection point
            subtract(grid, temp, perpendicular);    // vector from the spatial points to the cylinder axis
            distance = norm(perpendicular);         // distance between the points and vessel axis   
            normalize(perpendicular, distance);     // normalize the perpendicular vector
            // angle between the projected B0 and the vector from the projection point to the point
            phi_c = dot_product(perpendicular, B0_prj); // cos(phi)
            phi_2c2_1 = 2 * phi_c * phi_c - 1;      // cos(2*phi)
            // calculate the fieldmap from the vessel 
            distance = cyl_rad / distance;
            m_pFieldmap[p] += distance<=1 ? 2*M_PI * (1-m_Y)*m_dChi * (distance * distance) * phi_2c2_1 * theta_s2 : 2*M_PI * (1-m_Y)*m_dChi * (theta_c*theta_c - 1/3);
            m_pMask[p] = distance<=1.f ? m_pMask[p] : 1;
        }        
        bar.progress(c, m_num_cylinders);
    } 
    bar.finish();
    m_BVF = std::accumulate(m_pMask, m_pMask+res3, 0) * 100.0 / res3;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Fieldmaps generated successfully! BVF: " << m_BVF << "% Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";
}

void cylinder::print_info()
{
    shape::print_info();
    std::cout << "  Radius: " << m_radius << " um\n";
    std::cout << "  Number of cylinders: " << m_num_cylinders << "\n";
    std::cout << "  Orientation: " << m_orientation << " rad\n";
}

// -------------------------------------------------------------------------- //

bool cylinder::run()
{
    print_info();
    create_grid();
    generate_shapes();
    generate_mask_fieldmap();
    save();
    return true;
}