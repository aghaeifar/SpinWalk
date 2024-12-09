
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_base.cpp
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

#include "phantom_base.h"

namespace phantom
{
phantom_base::phantom_base()
{
    m_fov = 0;
    m_resolution = 0;
    m_seed = std::random_device{}();
    set_blood_parameters(0.273e-6 * 0.4, 0, 10.0);
    set_filename(); 
}

phantom_base::phantom_base(float fov_um, size_t resolution, float dChi, float Y, float BVF, int32_t seed, std::string filename)
:phantom_base()
{
    set_space(fov_um, resolution);
    set_blood_parameters(dChi, Y, BVF);
    set_filename(filename);
    if(seed >= 0)
        m_seed = seed;
}

phantom_base::~phantom_base()
{
}

void phantom_base::set_space(float fov_um, size_t resolution)
{
    this->m_fov = fov_um;
    this->m_resolution = resolution;
}

void phantom_base::set_blood_parameters(float dChi, float Y, float BVF)
{
    this->m_dChi = dChi;
    this->m_Y = Y;
    this->m_BVF = BVF;
    m_calc_fieldmap = Y >= 0;
}

void phantom_base::set_filename(std::string filename)
{
    this->m_filename = filename;
}

bool phantom_base::save() const
{
    BOOST_LOG_TRIVIAL(info) << "Saving the results..." << std::endl;
    std::filesystem::path parent_path = std::filesystem::absolute(m_filename).parent_path();
    if (std::filesystem::is_directory(parent_path) == false)
    {
        BOOST_LOG_TRIVIAL(warning) << "cannot find directory " << parent_path.string() << ". Trying to create it.";
        if(std::filesystem::create_directories(parent_path) == false)
        {
            BOOST_LOG_TRIVIAL(error)  << "cannot create directory " << parent_path.string();
            return false;
        }
    }

    HighFive::File file(m_filename, HighFive::File::Truncate);    
    // save fieldmap and mask
    std::vector<size_t> dims(3, m_resolution);
    if (m_Y >= 0)
    {
        HighFive::DataSet dataset_fieldmap = file.createDataSet<float>("fieldmap", HighFive::DataSpace(dims));
        dataset_fieldmap.write_raw((float *)m_fieldmap.data());
    }
    HighFive::DataSet dataset_mask = file.createDataSet<uint8_t>("mask", HighFive::DataSpace(dims));
    dataset_mask.write_raw((int8_t *)m_mask.data());
    // save fov
    float fov[3] = {m_fov*1e-6f, m_fov*1e-6f, m_fov*1e-6f}; // convert um to m
    std::vector<size_t> dims_fov(1, 3);
    HighFive::DataSet dataset_fov = file.createDataSet<float>("fov", HighFive::DataSpace(dims_fov));
    dataset_fov.write_raw(fov);
    // blood volume fraction
    std::vector<size_t> dims_1(1, 1);
    HighFive::DataSet dataset_BVF = file.createDataSet<float>("bvf", HighFive::DataSpace(dims_1));
    dataset_BVF.write_raw(&m_BVF);

    return true;
}

bool phantom_base::create_grid()
{
    BOOST_LOG_TRIVIAL(info) << "Creating grid..." << std::endl;
    if (m_fov == 0 || m_resolution == 0)
    {
        std::cerr << "Error: FOV or resolution is not set!" << std::endl;
        return false;
    }
    auto s = std::chrono::high_resolution_clock::now();
    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;
    m_grid.reserve(res3 * 3);
    // create a base grid
    std::vector<float> grid_base(res1, 0.f);
    double start = m_fov/res1/2.0;
    double end   = m_fov - m_fov/res1/2.0;
    double step  = (end - start) / (res1 - 1.0);
    for (int i = 0; i < res1; ++i) 
        grid_base[i] = start + i * step;
    // data is stored in row-major order in h5!!! 
    float *grid = (float *)m_grid.data();
    for (const auto& x : grid_base) 
        for (const auto& y : grid_base) 
            for (const auto& z : grid_base) 
            {
                grid[0] = x;
                grid[1] = y;
                grid[2] = z;
                grid += 3;
            }
    auto e = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Grid created successfully! " << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(e - s).count() << " s\n";
    return true;
}

// Overload << operator
std::ostream& operator<<(std::ostream& os, const phantom_base& obj) {
    os << "Phantom information:" << "\n" 
       << "  FOV: " << obj.m_fov << " \xC2\xB5m\n"
       << "  Resolution: " << obj.m_resolution << "\n"
       << "  Blood parameters: dChi=" << obj.m_dChi << ", Y=" << obj.m_Y <<  "\n"
       << "  Filename: " << obj.m_filename <<  "\n"
       << "  Seed: " << obj.m_seed <<  "\n";
    return os;
}



float dot_product(const float *a, const float *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

float norm(const float *a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

float norm_p2(const float *a)
{
    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]; // avoid sqrt
}

void cross_product(const float *a, const float *b, float *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

void normalize(float *a, float n)
{
    if (n < 0)
        n = norm(a);
    a[0] /= n;
    a[1] /= n;
    a[2] /= n;
}

void subtract(const float *a, const float *b, float *c)
{
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

void add(const float *a, const float *b, float *c)
{
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

void multiply(const float a, const float *b, float *c)
{
    c[0] = a * b[0];
    c[1] = a * b[1];
    c[2] = a * b[2];
}

void multiply(const float *a, const float *b, float *c)
{
    c[0] = a[0] * b[0];
    c[1] = a[1] * b[1];
    c[2] = a[2] * b[2];
}

void copy(const float *a, float *b)
{
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
}

void roty(float theta, const float *m0, float *m1)
{
    float deg2rad = 0.0174532925199433; // = M_PI/180 
    float s = sin(theta * deg2rad);
    float c = cos(theta * deg2rad);
    m1[0] =  c*m0[0] + s*m0[2];
    m1[1] =  m0[1];
    m1[2] = -s*m0[0] + c*m0[2];
}
}