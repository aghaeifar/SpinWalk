
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
#include <ranges>

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
    set_parameters(0.273e-6 * 0.4, 0, 10.0);
    set_filename(); 
}

phantom_base::phantom_base(float fov_um, size_t resolution, float dChi, float Y, float BVF, int32_t seed, std::string filename)
:phantom_base()
{
    set_space(fov_um, resolution);
    set_parameters(dChi, Y, BVF);
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

void phantom_base::set_parameters(float dChi, float Y, float volume_fraction)
{
    this->m_dChi = dChi;
    this->m_Y = Y;
    this->m_volume_fraction = volume_fraction;
    m_calc_fieldmap = Y >= 0;
}

void phantom_base::set_filename(std::string filename)
{
    this->m_filename = filename;
}

bool phantom_base::save() const
{
    float fov[3] = {m_fov*1e-6f, m_fov*1e-6f, m_fov*1e-6f}; // convert um to m
    std::vector<size_t> resolution(3, m_resolution);
    return save(fov, resolution);
}

bool phantom_base::save(float fov[3], std::vector<size_t> resolution) const
{
    BOOST_LOG_TRIVIAL(info) << "Saving the results...";
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
    if (m_Y >= 0)
    {
        HighFive::DataSet dataset_fieldmap = file.createDataSet<float>("fieldmap", HighFive::DataSpace(resolution));
        dataset_fieldmap.write_raw((float *)m_fieldmap.data());
    }
    HighFive::DataSet dataset_mask = file.createDataSet<uint8_t>("mask", HighFive::DataSpace(resolution));
    dataset_mask.write_raw((int8_t *)m_mask.data());
    // save fov
    std::vector<size_t> dims_fov(1, 3);
    HighFive::DataSet dataset_fov = file.createDataSet<float>("fov", HighFive::DataSpace(dims_fov));
    dataset_fov.write_raw(fov);
    // blood volume fraction
    std::vector<size_t> dims_1(1, 1);
    HighFive::DataSet dataset_BVF = file.createDataSet<float>("bvf", HighFive::DataSpace(dims_1));
    dataset_BVF.write_raw(&m_volume_fraction);

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
    const size_t resolution = m_resolution;
    const size_t resolutionSquared = resolution * resolution;
    const size_t resolutionCubed = resolution * resolutionSquared;
    // Reserve space for 3D coordinates (x, y, z for each point)
    m_grid.reserve(resolutionCubed * 3);
    // Create base grid using modern initialization
    std::vector<float> gridBase(resolution, 0.0f);
    // Calculate step parameters using const for immutable values
    const double start = m_fov / resolution / 2.0;
    const double end = m_fov - m_fov / resolution / 2.0;
    const double step = (end - start) / (resolution - 1.0);
    for (size_t i : std::views::iota(size_t(0), resolution)) 
        gridBase[i] = static_cast<float>(start + i * step);
    
    // data is stored in row-major order in h5!!! 
    float *grid = (float *)m_grid.data();
    for (const auto& x : gridBase) 
        for (const auto& y : gridBase) 
            for (const auto& z : gridBase) {
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

}