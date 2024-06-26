
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
#include "shape_base.cuh"


shape::shape()
{
    m_fov = 0;
    m_resolution = 0;
    m_random_seed = std::random_device{}();
    set_blood_parameters(0.273e-6 * 0.4, 0, 10.0);
    set_filename(); 
}

shape::shape(float fov_um, size_t resolution, float dChi, float Y, float BVF, bool is_seed_fixed, std::string filename)
:shape()
{
    set_space(fov_um, resolution);
    set_blood_parameters(dChi, Y, BVF);
    set_filename(filename);
    if(is_seed_fixed)
        m_random_seed = 0;
}

shape::~shape()
{
}

void shape::set_space(float fov_um, size_t resolution)
{
    this->m_fov = fov_um;
    this->m_resolution = resolution;
}

void shape::set_blood_parameters(float dChi, float Y, float BVF)
{
    this->m_dChi = dChi;
    this->m_Y = Y;
    this->m_BVF = BVF;
    m_calc_fieldmap = Y >= 0;
}

void shape::set_filename(std::string filename)
{
    this->m_filename = filename;
}

bool shape::save()
{
    std::cout << "Saving the results..." << std::endl;
    std::filesystem::path parent_path = std::filesystem::absolute(m_filename).parent_path();
    if (std::filesystem::is_directory(parent_path) == false)
    {
        std::cout << "cannot find directory " << parent_path.string() << ". Trying to create it.";
        if(std::filesystem::create_directories(parent_path) == false)
        {
            std::cerr  << "cannot create directory " << parent_path.string();
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

bool shape::create_grid()
{
    std::cout << "Creating grid..." << std::endl;
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
    std::cout << "Grid created successfully! " << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(e - s).count() << " s\n";
    return true;
}

void shape::print_info()
{
    std::cout << "Phantom information:" << "\n";
    std::cout << "  FOV: " << m_fov << " um\n";
    std::cout << "  Resolution: " << m_resolution << "\n";
    std::cout << "  Blood parameters: dChi=" << m_dChi << ", Y=" << m_Y <<  "\n";
    std::cout << "  Filename: " << m_filename <<  "\n";
}
