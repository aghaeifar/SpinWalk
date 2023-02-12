
#include <filesystem>
#include <stdexcept> // For standard exceptions
#include <cmath>

// boost includes
#include <boost/log/trivial.hpp> 

#include "pgse.h"
#include "ini.h"
#include "definitions.h"

namespace dMRI
{

std::string read(std::string config_file, std::string section, std::string key)
{
    BOOST_LOG_TRIVIAL(info) << "Reading config: " << config_file;
    if (std::filesystem::exists(config_file) == false)
        throw std::runtime_error( "Config-file does not exist: " + config_file);

    mINI::INIFile file(config_file);
    mINI::INIStructure ini;
    if(file.read(ini) == false)
        throw std::runtime_error("Failed to read config file: " + config_file );

    std::string value = "";
    if(ini["GENERAL"]["PARENT_CONFIG"].empty() == false) {
        std::filesystem::path parent_config(ini["GENERAL"]["PARENT_CONFIG"]);   
        if (parent_config.is_relative()) // if parent_config is relative, make it absolute
            parent_config = std::filesystem::absolute(config_file).parent_path() / parent_config;
        value = read(parent_config.string(), section, key);
    }
    if(ini[section][key].empty() == false)
        value = ini[section][key];

    return value;
}

void pgse::set_parameters(double b_value, uint32_t start_ms, uint32_t delta_ms, uint32_t DELTA_ms)
{
    this->b_value  = b_value;
    this->start_ms = start_ms;
    this->DELTA_ms = DELTA_ms;
    this->delta_ms = delta_ms; 
}

bool pgse::read_timestep(std::string config_file, int &timestep_us)
{
    std::string timestep_str;
    try{
        timestep_str = read(config_file, "SCAN_PARAMETERS", "TIME_STEP");
    } catch (const std::exception& e) {
        BOOST_LOG_TRIVIAL(error) << e.what();
        return false;
    }
    if (timestep_str.empty())
    {
        BOOST_LOG_TRIVIAL(error) << "TIME_STEP is not set! Create a section SCAN_PARAMETERS and set TIME_STEP";
        return false;
    }

    timestep_us = std::stoi(timestep_str);
    return true;
}

bool pgse::run(std::string config_file)
{
    int timestep_us;
    if(read_timestep(config_file, timestep_us) == false)
        return false;
    BOOST_LOG_TRIVIAL(info) << "TIME_STEP: " << timestep_us;

    if(DELTA_ms < delta_ms)
    {
        BOOST_LOG_TRIVIAL(error) << "\xCE\x94 must be greater than \xCE\xB4: " << DELTA_ms << " vs " << delta_ms;
        return false;
    }

    double d  = delta_ms * 1e-3; // gradient duration
    double D  = DELTA_ms * 1e-3; // distance between gradients
    double G2 = b_value  * 1e6 / (GAMMA * GAMMA * d * d * (D-d/3.0));
    double G  = sqrt(G2) * 1000.;  //  mT/m
    BOOST_LOG_TRIVIAL(info) << "Gradient amplitude: " << G << " mT/m";

    // normalize direction
    double norm = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    if (norm == 0)
    {
        BOOST_LOG_TRIVIAL(error) << "Direction vector is zero!";
        return false;
    }
    for(auto &d : dir)
        d = d / norm;

    mINI::INIFile file(config_file);
    mINI::INIStructure ini;
    file.read(ini); 

    uint32_t start_us = start_ms*1000; // time in us

    ini["SCAN_PARAMETERS"]["RF_FA[0]"] = "90";
    ini["SCAN_PARAMETERS"]["RF_FA[1]"] = "180";
    ini["SCAN_PARAMETERS"]["RF_PH[0]"] = "0";
    ini["SCAN_PARAMETERS"]["RF_PH[1]"] = "90";
    ini["SCAN_PARAMETERS"]["RF_T[0]"]  = "0";
    ini["SCAN_PARAMETERS"]["RF_T[1]"]  = std::to_string(start_us + delta_ms*1000 + (DELTA_ms-delta_ms)*1000/2);

    size_t i = 0, n_points = delta_ms * 1000 / timestep_us;
    std::string gradient_str = std::to_string(G * dir[0]) + " " + std::to_string(G * dir[1]) + " " + std::to_string(G * dir[2]);
    for (i=0; i < n_points; i++)
        ini["SCAN_PARAMETERS"]["GRADIENT_XYZ[" + std::to_string(i) + "]"] = gradient_str;
    for (; i < 2*n_points; i++)
        ini["SCAN_PARAMETERS"]["GRADIENT_XYZ[" + std::to_string(i) + "]"] = gradient_str;

    for (i=0; i < n_points; i++)
        ini["SCAN_PARAMETERS"]["GRADIENT_T[" + std::to_string(i) + "]"] = std::to_string(start_us + i*timestep_us);
    for (; i < 2*n_points; i++)
        ini["SCAN_PARAMETERS"]["GRADIENT_T[" + std::to_string(i) + "]"] = std::to_string(start_us + (DELTA_ms-delta_ms)*1000 + i*timestep_us);
    
    if(file.write(ini, true) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Failed to write config file: " << config_file;
        return false;
    }
    return true;
}

} // namespace dMRI

