
/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: miscellaneous.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#ifndef __CONFIG_READER_H__
#define __CONFIG_READER_H__

#include <map>
#include "ini.h"
#include "miscellaneous.h"


bool read_config(std::string  config_file, simulation_parameters& param, std::vector<float>& sample_length_scales, std::map<std::string, std::vector<std::string> >& filenames)
{
    if (std::filesystem::exists(config_file) == false)
    {
        std::cout << "File does not exist: " << config_file << std::endl;
        return false;
    }
    // first, create a file instance
    mINI::INIFile file(config_file);
    // next, create a structure that will hold data
    mINI::INIStructure ini;
    if (file.read(ini) == false)
    {
        std::cout << "Problem reading config file: " << config_file << std::endl;
        return false;
    }

    // reading section FILES
    if(ini.has("FILES"))
    {
        if (ini.get("FILES").has("FIELD_MAP[0]"))
        {
            filenames.at("fieldmap").clear();
            for (int i = 0; ini.get("FILES").has("FIELD_MAP[" + std::to_string(i) + "]"); i++)
            {
                filenames.at("fieldmap")
                    .push_back(ini.get("FILES").get("FIELD_MAP["
                        + std::to_string(i) + "]"));
                if (std::filesystem::exists(filenames.at("fieldmap").back())
                    == false)
                {
                    std::cout << "File does not exist: "
                        << filenames.at("fieldmap").back() << std::endl;
                    return false;
                }
            }
        }
        param.n_fieldmaps = filenames.at("fieldmap").size();

        if (ini.get("FILES").has("M0"))
        {
            filenames.at("M0").clear();
            filenames.at("M0").push_back(ini.get("FILES").get("M0"));
            if (std::filesystem::exists(filenames.at("M0").back()) == false)
            {
                std::cout << "File does not exist: " << filenames.at("M0").back()
                    << std::endl;
                return false;
            }
        }

        if (ini.get("FILES").has("OUTPUTS"))
        {
            filenames.at("output").clear();
            filenames.at("output").push_back(ini.get("FILES").get("OUTPUTS"));
        }
    }

    // reading section SCAN_PARAMETERS
    if(ini.has("SCAN_PARAMETERS"))
    {
        if(ini.get("SCAN_PARAMETERS").has("TR"))
            param.TR = std::stof(ini.get("SCAN_PARAMETERS").get("TR"));    
        if(ini.get("SCAN_PARAMETERS").has("DWELL_TIME"))
            param.dt = std::stof(ini.get("SCAN_PARAMETERS").get("DWELL_TIME"));
        if(ini.get("SCAN_PARAMETERS").has("DUMMY_SCAN"))
            param.n_dummy_scan  = std::stoi(ini.get("SCAN_PARAMETERS").get("DUMMY_SCAN"));
        if(ini.get("SCAN_PARAMETERS").has("FA"))
            param.FA = std::stof(ini.get("SCAN_PARAMETERS").get("FA")) * M_PI / 180.; // convert to radian
        
        param.TE = param.TR / 2.; // std::stof(ini.get("SCAN_PARAMETERS").get("TE"));
    }

    // reading section SIMULATION_PARAMETERS
    if(ini.has("SIMULATION_PARAMETERS"))
    {
        if(ini.get("SIMULATION_PARAMETERS").has("B0"))
            param.B0 = std::stof(ini.get("SIMULATION_PARAMETERS").get("B0"));
        if(ini.get("SIMULATION_PARAMETERS").has("SEED"))
            param.seed = std::stoi(ini.get("SIMULATION_PARAMETERS").get("SEED"));
        if(ini.get("SIMULATION_PARAMETERS").has("NUMBER_OF_SPINS"))
            param.n_spins = std::stof(ini.get("SIMULATION_PARAMETERS").get("NUMBER_OF_SPINS"));
        if(ini.get("SIMULATION_PARAMETERS").has("DIFFUSION_CONSTANT"))
            param.diffusion_const = std::stof(ini.get("SIMULATION_PARAMETERS").get("DIFFUSION_CONSTANT"));
        if(ini.get("SIMULATION_PARAMETERS").has("ENABLE_180_REFOCUSING"))
            param.enRefocusing180 = ini.get("SIMULATION_PARAMETERS").get("ENABLE_180_REFOCUSING").compare("0") != 0;
        if(ini.get("SIMULATION_PARAMETERS").has("SAMPLE_LENGTH_SCALES[0]"))
        {
            sample_length_scales.clear();
             for (int i = 0; ini.get("SIMULATION_PARAMETERS").has("SAMPLE_LENGTH_SCALES[" + std::to_string(i) + "]"); i++)
                sample_length_scales.push_back(std::stof(ini.get("SIMULATION_PARAMETERS").get("SAMPLE_LENGTH_SCALES[" + std::to_string(i) + "]")));
            param.n_sample_length_scales = sample_length_scales.size();
        }
    }

    // reading section TISSUE_PARAMETERS
    if(ini.has("TISSUE_PARAMETERS"))
    {
        if(ini.get("TISSUE_PARAMETERS").has("T1"))
            param.T1 = std::stof(ini.get("TISSUE_PARAMETERS").get("T1"));
        if(ini.get("TISSUE_PARAMETERS").has("T2"))
            param.T2 = std::stof(ini.get("TISSUE_PARAMETERS").get("T2"));
    }

    // reading section DEBUG 
    if(ini.has("DEBUG"))
    {
        if(ini.get("DEBUG").has("DUMP_INFO"))
            param.enDebug = ini.get("DEBUG").get("DUMP_INFO").compare("0") != 0;
        if(ini.get("DEBUG").has("SIMULATE_STEADYSTATE"))
            param.enSteadyStateSimulation  = ini.get("DEBUG").get("SIMULATE_STEADYSTATE").compare("0") != 0;
    }   

    param.prepare(); 
    return true;
}

#endif  // __CONFIG_READER_H__