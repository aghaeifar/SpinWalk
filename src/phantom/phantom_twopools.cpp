
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_twopools.cpp
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
#include "phantom_twopools.h"

// -------------------------------------------------------------------------- //
namespace phantom
{


twopools::twopools()
{

}

twopools::twopools(float fov_um, size_t resolution, std::string filename)
: phantom_base(fov_um, resolution, 0, -1, 0, 0, filename)
{

}

twopools::~twopools()
{
}


bool twopools::generate_mask_fieldmap()
{
    BOOST_LOG_TRIVIAL(info) << "Generating twopools..." << std::endl;
    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;
    int32_t x_min, x_max, y_min, y_max, z_min, z_max;

    auto start = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Allocating memory...";
    m_fieldmap.resize(0, 0.f);
    m_mask.resize(res3, 0);
    BOOST_LOG_TRIVIAL(info) << "Done!\n";

    std::fill(m_mask.begin(), m_mask.begin() + m_mask.size()/2, 1);

    m_volume_fraction = std::accumulate(m_mask.begin(), m_mask.end(), 0) * 100.0 / m_mask.size();
    BOOST_LOG_TRIVIAL(info) << "Actual Volume Fraction = " << m_volume_fraction << "% ...\n";   
    auto end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Twopools phantom generated successfully! Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";
    return true;
}


std::ostream& operator<<(std::ostream& os, const twopools& obj)
{
    os << static_cast<const phantom_base&>(obj);
    return os;
}

// -------------------------------------------------------------------------- //

bool twopools::run(bool write_to_disk)
{
    BOOST_LOG_TRIVIAL(info) << *this << std::endl;
    if(create_grid() == false)
        return false;
    if(generate_mask_fieldmap() == false)
        return false;
    if(write_to_disk)
        if(save() == false)
            return false;
    return true;
}

}