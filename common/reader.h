
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

namespace reader
{

bool read_config(std::string config_filename, simulation_parameters& param, std::vector<float>& sample_length_scales, std::map<std::string, std::vector<std::string> >& filenames)
{
    if (std::filesystem::exists(config_filename) == false)
    {
        std::cout << ERR_MSG << "file does not exist: " << config_filename << std::endl;
        return false;
    }
    // first, create a file instance
    mINI::INIFile file(config_filename);
    // next, create a structure that will hold data
    mINI::INIStructure ini;
    if (file.read(ini) == false)
    {
        std::cout << ERR_MSG << "problem reading config file: " << config_filename << std::endl;
        return false;
    }

    if(ini.has("parent"))
    {
        if(ini.get("parent").has("parent_config"))
        {
            std::string parent_config = ini.get("parent").get("parent_config");            
            parent_config = std::filesystem::absolute(parent_config);
            std::cout << "Reading parent config: " << parent_config << std::endl;
            if (read_config(parent_config, param, sample_length_scales, filenames) == false)
                return false;
        }
    }

    // reading section FILES
    if(ini.has("files"))
    {   // only read fieldmap to count number of fieldmaps and set n_fieldmaps
        if(ini.get("files").has("fieldmap[0]"))
        {
            filenames.at("fieldmap").clear();
            for (uint16_t i = 0; ini.get("files").has("fieldmap[" + std::to_string(i) + "]"); i++)                
                filenames.at("fieldmap").push_back(ini.get("files").get("fieldmap[" + std::to_string(i) + "]"));
        }
        param.n_fieldmaps = filenames.at("fieldmap").size();

        // read all other files 
        for (std::map<std::string, std::vector<std::string>>::iterator it=filenames.begin(); it!=filenames.end(); ++it)
        {
            if( it->first.compare("fieldmap") == 0)
                continue; // fieldmap is already read
            it->second.clear();
            for (uint16_t i = 0; i<param.n_fieldmaps; i++)
                if(ini.get("files").has(it->first + "[" + std::to_string(i) + "]"))
                    it->second.push_back(ini.get("files").get(it->first + "[" + std::to_string(i) + "]"));
                else
                    it->second.push_back("");
        }

        // if m1 is empty, create filename for m1. filename is the same as fieldmap but with prefix config_filename
        for(uint16_t i=0; i<param.n_fieldmaps; i++)
            if(filenames.at("m1")[i].empty())
            {             
                std::string f_config  = std::filesystem::path(config_filename).filename().string();
                std::string f_field   = std::filesystem::path(filenames.at("fieldmap")[i]).filename().string();
                std::string parent    = std::filesystem::path(filenames.at("fieldmap")[i]).parent_path().string();
                std::string ext       = std::filesystem::path(filenames.at("fieldmap")[i]).extension().string();
                filenames.at("m1")[i] = parent + "/" + f_config + "_" + f_field + ext;
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
            param.FA = std::stof(ini.get("SCAN_PARAMETERS").get("FA")) ; // convert to radian
        if(ini.get("SCAN_PARAMETERS").has("APPLY_FA/2"))
            param.enApplyFA2 = ini.get("SCAN_PARAMETERS").get("APPLY_FA/2").compare("0") != 0;
        if(ini.get("SCAN_PARAMETERS").has("PHASE_CYCLING"))
            param.phase_cycling = std::stof(ini.get("SCAN_PARAMETERS").get("PHASE_CYCLING")) ; // convert to radian

        if(ini.get("SCAN_PARAMETERS").has("ENABLE_REFOCUSING"))
            param.enRefocusing = ini.get("SIMULATION_PARAMETERS").get("ENABLE_REFOCUSING").compare("0") != 0;
        
        uint16_t i = 0;
        for(i=0; i<MAX_TE && ini.get("SCAN_PARAMETERS").has("TE[" + std::to_string(i) + "]"); i++)
        {
            param.TE[i] = std::stof(ini.get("SCAN_PARAMETERS").get("TE[" + std::to_string(i) + "]")) / param.dt;
            if (param.TE[i] == 0)
                break;
        }
        param.n_TE = i;
        if(param.n_TE == 0)
        {
            std::cout << ERR_MSG << "TE is not defined" << std::endl;
            return false;
        }

        for(i=0; i<MAX_SE && ini.get("SCAN_PARAMETERS").has("RF_SE[" + std::to_string(i) + "]"); i++)
            param.RF_SE[i] = std::stof(ini.get("SCAN_PARAMETERS").get("RF_SE[" + std::to_string(i) + "]")) ; // convert to radian
        int ii = i;

        for(i=0; i<MAX_SE && ini.get("SCAN_PARAMETERS").has("RF_SE_PHS[" + std::to_string(i) + "]"); i++)
            param.RF_SE_PHS[i] = std::stof(ini.get("SCAN_PARAMETERS").get("RF_SE_PHS[" + std::to_string(i) + "]")) ; // convert to radian

        if(ii != i)
        {
            std::cout << ERR_MSG << "RF_SE and RF_SE_PHS must have the same number of elements" << std::endl;
            return false;
        }

        for(i=0; i<MAX_SE && ini.get("SCAN_PARAMETERS").has("T_SE[" + std::to_string(i) + "]"); i++)
        {
            param.T_SE[i] = std::stof(ini.get("SCAN_PARAMETERS").get("T_SE[" + std::to_string(i) + "]")) / param.dt;
            if (param.T_SE[i] == 0)
                break;
        }
        if(ii != i)
        {
            std::cout << ERR_MSG << "Not enought or too many T_SE is defined (expected " << ii << " elements but found " << i << ")" << std::endl;
            return false;
        }

        param.n_SE = i;
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

bool read_header(std::string filename, input_header &hdr_in)
{
    if(std::filesystem::exists(filename) == false)
    {
        std::cout << ERR_MSG << "file does not exist: " << filename << std::endl;
        return false;
    }

    std::ifstream in_field(filename, std::ios::in | std::ios::binary);
    if (!in_field.is_open()) 
    {
        std::cout << ERR_MSG << "problem reading header of " << filename << std::endl;
        return false;
    }

    in_field.read((char*)&hdr_in, sizeof(input_header));
    return true;
}

bool read_fieldmap(std::string fieldmap_filename, std::vector<float> &fieldmap, std::vector<char> &mask, simulation_parameters& param)
{
    if(std::filesystem::exists(fieldmap_filename) == false)
    {
        std::cout << ERR_MSG << "fieldmap file does not exist: " << fieldmap_filename << std::endl;
        return false;
    }

    input_header hdr_in;
    std::cout << "Opening fieldmap " << fieldmap_filename << std::endl;
    std::ifstream in_field(fieldmap_filename, std::ios::in | std::ios::binary);
    if (!in_field.is_open()) 
    {
        std::cout << ERR_MSG << "problem opening file " << fieldmap_filename << std::endl;
        return false;
    }

    in_field.read((char*)&hdr_in, sizeof(input_header));
    hdr_in.print();
    std::copy(hdr_in.fieldmap_size, hdr_in.fieldmap_size + 3, param.fieldmap_size);
    std::copy(hdr_in.sample_length, hdr_in.sample_length + 3, param.sample_length);
    param.matrix_length = param.fieldmap_size[0] * param.fieldmap_size[1] * param.fieldmap_size[2];
    if (fieldmap.size() != param.matrix_length)
    {
        std::cout << "Fieldmap size changed. Re-allocating memory..." << std::endl;
        std::cout << "Old size: " << fieldmap.size() << std::endl;
        std::cout << "New size: " << param.matrix_length << std::endl;
        std::cout << "New length (um): " << param.sample_length[0] * 1e6 << " " << param.sample_length[1] * 1e6 << " " << param.sample_length[2] * 1e6 << std::endl;
        fieldmap.resize(param.matrix_length);
        mask.resize(param.matrix_length);
    }

    std::cout << "Reading...fieldmap...";
    in_field.read((char*)fieldmap.data(), sizeof(float) * param.matrix_length); std::cout << "done...mask...";
    in_field.read((char*)mask.data(), sizeof(bool) * param.matrix_length); std::cout << "done." << std::endl;
    in_field.close();
    return true;
}

template <typename T>
bool read_file(std::string filename, std::vector<T> &storage)
{
    if(std::filesystem::exists(filename) == false)
    {
        std::cout << ERR_MSG << "file does not exist: " << filename << std::endl;
        return false;
    }

    std::cout << "Reading " << filename << std::endl;
    std::ifstream in_field(filename, std::ios::in | std::ios::binary);
    if (!in_field.is_open()) 
    {
        std::cout << ERR_MSG << "error opening file " << filename << std::endl;
        return false;
    }
    
    in_field.read((char*)storage.data(), sizeof(storage[0]) * storage.size());
    in_field.close();

    return true;
}

}
#endif  // __CONFIG_READER_H__