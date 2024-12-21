
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: arbitrary_gradient.h
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
#include "phantom_cylinder.h"


namespace phantom
{

bool check_cylinder_overlap(const std::vector<std::vector<float>>& m_cylinder_points, const std::vector<float>& m_cylinder_radii, const float* cyl_pnt, float& radius, bool is_random_radius) {
    int c = 0;
    bool stop = false;

    #pragma omp parallel for shared(stop, radius)
    for (c = 0; c < m_cylinder_radii.size(); c++) {
        if (stop) continue; // Skip computation if condition already met.

        float p2p1[3];
        subtract(cyl_pnt, m_cylinder_points[c].data(), p2p1);
        float distance = sqrtf(p2p1[0]*p2p1[0] + p2p1[1]*p2p1[1]);;
        // Check the conditions
        if (distance <= m_cylinder_radii[c] || distance <= radius) {
            #pragma omp critical
            {
                stop = true; // Signal all threads to stop further checks.
            }
        } else if (distance < m_cylinder_radii[c] + radius) {
            if (!is_random_radius) {
                #pragma omp critical
                {
                    stop = true; // Signal all threads to stop further checks.
                }
            } else {
                #pragma omp critical
                {
                    if (!stop) 
                        radius = distance - m_cylinder_radii[c]; // Adjust radius.
                }
            }
        }
    }
    return stop; // Return whether an overlap condition was found.
}


cylinder::cylinder()
{
    m_radius = 0;
    m_orientation = 0;
}

cylinder::cylinder(float fov_um, size_t resolution, float dChi, float Y, float radius_um, float volume_fraction, float orientation, int32_t seed, std::string filename)
: phantom_base(fov_um, resolution, dChi, Y, volume_fraction, seed, filename)
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
    roty(orientation, B0_orig, B0);
}

bool cylinder::generate_shapes()
{
    if(2*m_radius>=m_fov)
    {
        BOOST_LOG_TRIVIAL(error) << "Error: The radius of the cylinder is too large for the given FOV!\n";
        return false;
    }
    BOOST_LOG_TRIVIAL(info) << "Generating coordinates...for target BVF = " << m_volume_fraction << "% ...\n"; 
    bool is_random_radius = m_radius < 0;
    float max_radius    = m_radius>0 ? m_radius:-m_radius;
    m_cylinder_points.clear();
    m_cylinder_radii.clear();
    float cyl_pnt[3], cyl_rad ;
    // srandom engine
    std::mt19937 gen(m_seed); // Mersenne Twister generator
    std::uniform_real_distribution<float> dist(0.f, 1.f); 
      
    float distance, vol_cyl = 0, vol_cyl_total = 0, vol_tol = m_fov*m_fov*m_fov;
    int32_t progress = 0;
    auto bar = barkeep::ProgressBar(&progress, {.total = 100, .message = "Generating (1/2)", .style = barkeep::ProgressBarStyle::Rich,});
    auto start = std::chrono::high_resolution_clock::now();
    while(progress < 100)
    {
        cyl_rad = is_random_radius ? dist(gen) * max_radius : max_radius;
        for (size_t i = 0; i < 3; i++) // generate a random point for a cylinder which fit in the FOV
            cyl_pnt[i] = dist(gen) * (m_fov+2*cyl_rad) - cyl_rad;   
        
        // check if cylinder coordinate is ok
        if (check_cylinder_overlap(m_cylinder_points, m_cylinder_radii, cyl_pnt, cyl_rad, is_random_radius)) 
            continue;   

        vol_cyl = calculate_volume(cyl_pnt, cyl_rad);
        // if the total volume of the cylinders is more than the target BVF or the cylinder is outside of volume, skip this cylinder
        if (100*(vol_cyl + vol_cyl_total) / vol_tol > 1.02*m_volume_fraction || vol_cyl < 0)
            continue;

        vol_cyl_total += vol_cyl;
        progress = 100*(100.*vol_cyl_total/vol_tol/m_volume_fraction);
        m_cylinder_points.push_back({cyl_pnt[0], cyl_pnt[1], cyl_pnt[2]});
        m_cylinder_radii.push_back(cyl_rad);     
    }
    bar->done();
    auto end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << m_cylinder_radii.size() << " coordinates generated successfully! Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << std::endl;
    return true;
}


float cylinder::calculate_volume(float *cyl_pnt, float cyl_rad)
{
    bool intersect = false;
    // check if the cylinder is completely outside the volume
    for (size_t i = 0; i < 2; i++)
        if(cyl_pnt[i]+cyl_rad < 0 || cyl_pnt[i]-cyl_rad > m_fov)
            return -1.f;
    // check if the cylinder is completely inside the volume
    for (size_t i = 0; i < 2; i++) // only x and y directions
        if (cyl_pnt[i] < cyl_rad - 1.5 || cyl_pnt[i] > m_fov - cyl_rad + 1.5) // 1.5 is a small margin because of a possible larger volume after discretization
            intersect = true;

    if (intersect == false)
        return M_PI * cyl_rad*cyl_rad * m_fov;    

    // find the bounding box of the cylinder
    float v_size = m_fov / m_resolution;
    float cyl_rad2  = cyl_rad*cyl_rad;
    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    int32_t cyl_pnt_vox[3] = {int32_t(cyl_pnt[0]/v_size), int32_t(cyl_pnt[1]/v_size), int32_t(cyl_pnt[2]/v_size)};
    int32_t cyl_rad_vox = std::ceil(cyl_rad / m_fov * m_resolution)+1; 
    int32_t z_min = 0;
    int32_t z_max = m_resolution;
    int32_t x_min = std::max(0, cyl_pnt_vox[0] - cyl_rad_vox);
    int32_t x_max = std::min((int32_t)m_resolution, cyl_pnt_vox[0] + cyl_rad_vox + 2);
    int32_t y_min = std::max(0, cyl_pnt_vox[1] - cyl_rad_vox);
    int32_t y_max = std::min((int32_t)m_resolution, cyl_pnt_vox[1] + cyl_rad_vox + 2);
    int32_t counter = 0;
    #pragma omp parallel for
    for(int32_t pz=z_min; pz<z_max; pz++)
    for(int32_t py=y_min; py<y_max; py++)
    for(int32_t px=x_min; px<x_max; px++)
    {
        size_t p = px*res2 + py*res1 + pz;
        float *grid = &m_grid[3*p];
        float p2p1[3];
        // distance between the points and vessel axis and vector from the projection point to the point
        subtract(grid, cyl_pnt, p2p1);  // vector from the spatial points to the cylinder point
        float distance2 = p2p1[0]*p2p1[0] + p2p1[1]*p2p1[1];   // distance^2 between the points and vessel axis. this is distance in the plane perpendicular to the vessel axis
        if (distance2 <= cyl_rad2)
        {
            #pragma omp atomic
            ++counter;
        }
    }
    return counter * v_size*v_size*v_size ;
}


bool cylinder::generate_mask_fieldmap()
{   
    // set the cylinders orientation if they are parallel
    BOOST_LOG_TRIVIAL(info) << "Generating cylinders..." << std::endl;    
    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;
    int32_t x_min, x_max, y_min, y_max, z_min, z_max;
    int32_t cyl_rad_vox, cyl_pnt_vox[3];

    BOOST_LOG_TRIVIAL(info) <<"B0 direction: ["<<B0[0]<<", "<<B0[1]<<", "<<B0[2]<<"]\n";
    BOOST_LOG_TRIVIAL(info) <<"Allocating memory...";
    m_fieldmap.resize(m_calc_fieldmap ? res3:0, 0.f);
    m_mask.resize(res3, 0);
    BOOST_LOG_TRIVIAL(info) << "Done!\n";
    float v_size = m_fov / m_resolution;

    float cyl_dir[3] = {0.0, 0.0, 1.0};
    float theta_c, theta_c2, theta_s2, B0_prj[3] = {B0[0], B0[1], 0.0};     // project B0 from the projection point to the point to the plane perpendicular to the vessel axis  
    normalize(B0_prj);

    theta_c  = cos(m_orientation * M_PI / 180); // cos(theta), angle between axis of vessel and B0 (in radian)
    theta_c2 = theta_c * theta_c;
    theta_s2 = 1. - theta_c2; // sin^2(theta)

    size_t c = 0;
    auto bar = barkeep::ProgressBar(&c, {.total = m_cylinder_radii.size(), .message = "Generating (2/2)", .style = barkeep::ProgressBarStyle::Rich,});
    auto start = std::chrono::high_resolution_clock::now();
    for (c = 0; c < m_cylinder_radii.size(); c++)
    {
        float *cyl_pnt  = m_cylinder_points[c].data();
        float cyl_rad   = m_cylinder_radii[c];
        float cyl_rad2  = cyl_rad*cyl_rad;
        cyl_rad_vox = std::ceil(cyl_rad / v_size)+1;
        cyl_pnt_vox[0] = int32_t(cyl_pnt[0]/v_size); 
        cyl_pnt_vox[1] = int32_t(cyl_pnt[1]/v_size);
        cyl_pnt_vox[2] = int32_t(cyl_pnt[2]/v_size);

        // find the bounding box of the cylinder
        z_min = 0;
        z_max = m_resolution;
        if (m_calc_fieldmap)
        {
            x_min = std::max(0, cyl_pnt_vox[0] - cyl_rad_vox*8);
            x_max = std::min((int32_t)m_resolution, cyl_pnt_vox[0] + cyl_rad_vox*10);
            y_min = std::max(0, cyl_pnt_vox[1] - cyl_rad_vox*8);
            y_max = std::min((int32_t)m_resolution, cyl_pnt_vox[1] + cyl_rad_vox*10);
        } else 
        {
            x_min = std::max(0, cyl_pnt_vox[0] - cyl_rad_vox);
            x_max = std::min((int32_t)m_resolution, cyl_pnt_vox[0] + cyl_rad_vox + 2);
            y_min = std::max(0, cyl_pnt_vox[1] - cyl_rad_vox);
            y_max = std::min((int32_t)m_resolution, cyl_pnt_vox[1] + cyl_rad_vox + 2);
        }

        #pragma omp parallel for
        for(int32_t pz=z_min; pz<z_max; pz++)
        for(int32_t py=y_min; py<y_max; py++)
        for(int32_t px=x_min; px<x_max; px++) 
        {
            size_t p = px*res2 + py*res1 + pz;
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
                // reference: https://doi.org/10.1016/j.neuroimage.2017.09.015
                if (distance2 > cyl_rad2)  // outside the cylinder              
                    m_fieldmap[p] += 2*M_PI * (1-m_Y)*m_dChi * (cyl_rad2 / distance2) * phi_2c2_1 * theta_s2;   
                else // inside the cylinder                
                    m_fieldmap[p] += 2*M_PI * (1-m_Y)*m_dChi * (theta_c2 - 1.0/3.0);
            }                 
        }        
    } 
    bar->done();    
    m_volume_fraction = std::accumulate(m_mask.begin(), m_mask.end(), 0) * 100.0 / m_mask.size();
    BOOST_LOG_TRIVIAL(info) << "Actual Volume Fraction = " << m_volume_fraction << "% ...\n";   
    auto end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Cylinders generated successfully! " << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";
    return true;
}

std::ostream& operator<<(std::ostream& os, const cylinder& obj)
{
    os << static_cast<const phantom_base&>(obj) 
       << "  Radius: " << obj.m_radius << " \xC2\xB5m\n"
       << "  Volume Fraction: " << obj.m_volume_fraction << "\n"
       << "  Orientation: " << obj.m_orientation << " rad\n";
    return os;
}

// -------------------------------------------------------------------------- //

bool cylinder::run()
{
    BOOST_LOG_TRIVIAL(info) << *this << std::endl;
    if(create_grid() == false)
        return false;
    if(generate_shapes() == false)
        return false;
    if(generate_mask_fieldmap() == false)
        return false;
    if(save() == false)
        return false;
    return true;
}

} // namespace phantom