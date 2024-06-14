
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
}

cylinder::cylinder(float fov_um, size_t resolution, float dChi, float Y, float radius_um, float BVF, float orientation, std::string filename)
: shape(fov_um, resolution, dChi, Y, BVF, filename)
{
    set_cylinder_parameters(radius_um, orientation);
}

cylinder::~cylinder()
{

}

void cylinder::set_cylinder_parameters(float radius, float orientation)
{
    m_radius = radius;
    m_orientation = orientation;
    float B0_orig[3] = {0.f, 0.f, 1.f};
    yrot(orientation, B0_orig, B0);
}

void cylinder::generate_shapes()
{
    std::cout << "Generating coordinates...for BVF = " << m_BVF << "% ...\n";    
    bool is_random_radius = m_radius < 0;
    float max_radius    = m_radius>0 ? m_radius:-m_radius;
    m_cylinder_points.clear();
    m_cylinder_radii.clear();
    float cyl_pnt[3], radius ;

    // srandom engine
    std::random_device rd; // Seed
    std::mt19937 gen(rd()); // Mersenne Twister generator
    std::uniform_real_distribution<> dist(0.f, 1.f); 
      
    float distance, vol_cyl = 0, vol_tol = m_fov*m_fov*m_fov;
    int32_t progress = 0, prg_l = 0;
    auto start = std::chrono::high_resolution_clock::now();
    tqdm bar;
    while(progress < 100)
    {
        radius = is_random_radius ? dist(gen) * max_radius : max_radius;
        for (size_t i = 0; i < 3; i++) // generate a random point for a sphere which fit in the FOV
            cyl_pnt[i] = dist(gen) * m_fov;   
        
        // check if sphere coordinate is ok
        size_t c;
        for (c=0; c<m_cylinder_radii.size(); c++)
        {   
            float p2p1[3];
            subtract(cyl_pnt, &m_cylinder_points[3*c], p2p1);
            distance = sqrtf(p2p1[0]*p2p1[0] + p2p1[1]*p2p1[1]);
            // if the sphere is inside another sphere, generate a new sphere
            if(distance <= m_cylinder_radii[c] ||  distance <= radius)
                break;
            // adjust the radius of the sphere to avoid overlap
            if (distance < m_cylinder_radii[c] + radius)
            {
                if (!is_random_radius)
                    break;            
                radius = distance - m_cylinder_radii[c];
            }
        }
        if (c < m_cylinder_radii.size())
            continue;

        vol_cyl += M_PI * radius*radius * m_fov; // volume of the cylinder
        m_cylinder_points.insert(m_cylinder_points.end(), cyl_pnt, cyl_pnt+3);
        m_cylinder_radii.push_back(radius);  
        progress = 0.95 * 100*(100.*vol_cyl/vol_tol/m_BVF); // 0.95 is a factor to compensate for cylinders in the boundary   
        if (progress >= prg_l)
        { 
            bar.progress(progress, 100);
            prg_l += 5;
        }
    }
    bar.finish();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << m_cylinder_radii.size() << " coordinates generated successfully! Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << std::endl;
}

void cylinder::generate_mask_fieldmap()
{   
    // set the cylinders orientation if they are parallel

    std::cout << "Generating cylinders..." << std::endl;    
    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;
    int32_t x_min, x_max, y_min, y_max, z_min, z_max;
    int32_t cyl_rad_vox, cyl_pnt_vox[3];

    std::cout<<"B0 direction: ["<<B0[0]<<", "<<B0[1]<<", "<<B0[2]<<"]\n";
    std::cout<<"Allocating memory..."<<std::endl;
    m_fieldmap.resize(m_calc_fieldmap ? res3:0, 0.f);
    m_mask.resize(res3, 0);
    float v_size = m_fov / m_resolution;

    float cyl_dir[3] = {0.0, 0.0, 1.0};
    float theta_c, theta_c2, theta_s2, B0_prj[3] = {B0[0], B0[1], 0.0};     // project B0 from the projection point to the point to the plane perpendicular to the vessel axis  
    normalize(B0_prj);

    theta_c  = cos(m_orientation * M_PI / 180); // cos(theta), angle between axis of vessel and B0 (in radian)
    theta_c2 = theta_c * theta_c;
    theta_s2 = 1. - theta_c2; // sin^2(theta)
    
    tqdm bar;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t c = 0; c < m_cylinder_radii.size(); c++)
    {
        float *cyl_pnt = &m_cylinder_points[3*c];
        float cyl_rad   = m_cylinder_radii[c];
        float cyl_rad2  = cyl_rad*cyl_rad;
        cyl_rad_vox = std::ceil(cyl_rad / v_size)+1;
        cyl_pnt_vox[0] = int32_t(cyl_pnt[0]/v_size); 
        cyl_pnt_vox[1] = int32_t(cyl_pnt[1]/v_size);
        cyl_pnt_vox[2] = int32_t(cyl_pnt[2]/v_size);

        // find the bounding box of the sphere
        z_min = 0;
        z_max = m_resolution;
        if (m_calc_fieldmap)
        {
            x_min = std::max(0, cyl_pnt_vox[0] - cyl_rad_vox*8);
            x_max = std::min((int32_t)m_resolution, cyl_pnt_vox[0] + cyl_rad_vox*8);
            y_min = std::max(0, cyl_pnt_vox[1] - cyl_rad_vox*8);
            y_max = std::min((int32_t)m_resolution, cyl_pnt_vox[1] + cyl_rad_vox*8);
        } else 
        {
            x_min = std::max(0, cyl_pnt_vox[0] - cyl_rad_vox);
            x_max = std::min((int32_t)m_resolution, cyl_pnt_vox[0] + cyl_rad_vox + 2);
            y_min = std::max(0, cyl_pnt_vox[1] - cyl_rad_vox);
            y_max = std::min((int32_t)m_resolution, cyl_pnt_vox[1] + cyl_rad_vox + 2);
        }

        #pragma omp parallel for
        for(int32_t px=x_min; px<x_max; px++)
        for(int32_t py=y_min; py<y_max; py++)
        for(int32_t pz=z_min; pz<z_max; pz++)
        {
            int p = px*res2 + py*res1 + pz;
            float *grid = &m_grid[3*p];
            float p2p1[3], temp[3], perpendicular[3], distance2, phi_c, phi_2c2_1 ;
            // distance between the points and vessel axis and vector from the projection point to the point
            subtract(grid, cyl_pnt, p2p1);  // vector from the spatial points to the cylinder point
            distance2 = p2p1[0]*p2p1[0] + p2p1[1]*p2p1[1];   // distance^2 between the points and vessel axis. this is distance in the plane perpendicular to the vessel axis
            if (distance2 <= cyl_rad2)
                m_mask[p] = 1;

            if (m_calc_fieldmap)
            {
                multiply(dot_product(cyl_dir, p2p1), cyl_dir, temp); // project vector temp onto the cylinder direction vector
                add(temp, cyl_pnt, temp);               // projection point
                subtract(grid, temp, perpendicular);    // vector from the spatial points to the cylinder axis
                // angle between the projected B0 and the vector from the projection point to the point
                phi_c = dot_product(perpendicular, B0_prj) / norm(perpendicular); // cos(phi)
                phi_2c2_1 = 2 * phi_c * phi_c - 1;      // cos(2*phi)
                // calculate the fieldmap from the vessel 
                if (distance2 > cyl_rad2)  // outside the cylinder              
                    m_fieldmap[p] += 2*M_PI * (1-m_Y)*m_dChi * (cyl_rad2 / distance2) * phi_2c2_1 * theta_s2;   
                else // inside the cylinder                
                    m_fieldmap[p] += 2*M_PI * (1-m_Y)*m_dChi * (theta_c2 - 1/3);
            }                 
        }        
        bar.progress(c, m_cylinder_radii.size());
    } 
    bar.finish();
    m_BVF = std::accumulate(m_mask.begin(), m_mask.end(), 0) * 100.0 / res3;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Cylinders generated successfully! BVF: " << m_BVF << "% Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";
}

void cylinder::print_info()
{
    shape::print_info();
    std::cout << "  Radius: " << m_radius << " um\n";
    std::cout << "  BVF: " << m_BVF << "\n";
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