
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: file_utils.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#ifndef __FILE_UTILS_H__
#define __FILE_UTILS_H__

#include <filesystem>
#include <vector>
#include <map>
#include <algorithm> 
#include "miscellaneous.h"
#include "ini.h"

namespace file_utils
{

template <typename T>
int sort_remove_duplicates(T *array, int n)
{
    std::vector<T> v(array, array+n);
    std::sort(v.begin(), v.end()); 
    v.erase( std::unique( v.begin(), v.end() ), v.end() );
    std::copy(v.begin(), v.end(), array);
    return v.size();
}

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
            parent_config = std::filesystem::absolute(parent_config).string();
            std::cout << "Reading parent config: " << parent_config << std::endl;
            if (read_config(parent_config, param, sample_length_scales, filenames) == false)
                return false;
        }
    }

    // ============== reading section FILES ==============
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
            it->second.assign(param.n_fieldmaps, ""); // this clear data from parent config files, is right?
            for (uint16_t i = 0; i<param.n_fieldmaps && ini.get("files").has(it->first + "[" + std::to_string(i) + "]"); i++)
                it->second[i] = ini.get("files").get(it->first + "[" + std::to_string(i) + "]");
        }

        // This is a mandatory output. If m1 is empty, create filename for m1. filename is the same as fieldmap but with prefix config_filename
        for(uint16_t i=0; i<param.n_fieldmaps; i++)
            if(filenames.at("m1")[i].empty())
            {             
                std::string f_config  = std::filesystem::path(config_filename).filename().string();
                std::string f_field   = std::filesystem::path(filenames.at("fieldmap")[i]).filename().string();
                std::string parent    = std::filesystem::path(filenames.at("fieldmap")[i]).parent_path().string();
                std::string ext       = std::filesystem::path(filenames.at("fieldmap")[i]).extension().string();
                filenames.at("m1")[i] = parent + "/" + f_config + "_" + f_field + "_m1" + ext;
            }    
    }


    // ============== reading section SCAN_PARAMETERS ==============
    if(ini.has("SCAN_PARAMETERS"))
    {
        if(ini.get("SCAN_PARAMETERS").has("TR"))
            param.TR = std::stof(ini.get("SCAN_PARAMETERS").get("TR"));    
        if(ini.get("SCAN_PARAMETERS").has("DWELL_TIME"))
            param.dt = std::stof(ini.get("SCAN_PARAMETERS").get("DWELL_TIME"));             

        uint16_t i = 0;
        // ---------------- Echo times ----------------       
        param.n_TE = 0;
        for(i=0; i<MAX_TE && ini.get("SCAN_PARAMETERS").has("TE[" + std::to_string(i) + "]"); i++)        
            param.TE[i] = std::stof(ini.get("SCAN_PARAMETERS").get("TE[" + std::to_string(i) + "]")) / param.dt;

        // check TE conditions
        if (std::is_sorted(param.TE, param.TE + i) == false || 
            std::adjacent_find(param.TE, param.TE + i) != param.TE + i || 
            param.TE[0] < 0 || 
            (param.n_TE = i) == 0)
        {
            std::cout << ERR_MSG << "TE must be in ascending order and must not have duplicates or negative values" << std::endl;
            return false;
        }        

        // ---------------- RF pulses (start times, Flip angles, phases and ) ----------------
        // RF start times
        for(i=0; i<MAX_RF && ini.get("SCAN_PARAMETERS").has("RF_ST[" + std::to_string(i) + "]"); i++)        
            param.RF_ST[i] = std::stof(ini.get("SCAN_PARAMETERS").get("RF_ST[" + std::to_string(i) + "]")) / param.dt;
        
        // check RF start time conditions
        if (std::is_sorted(param.RF_ST, param.RF_ST + i) == false || 
            std::adjacent_find(param.RF_ST, param.RF_ST + i) != param.RF_ST + i || 
            param.RF_ST[0] != 0 || 
            (param.n_RF = i) == 0)
        {
            std::cout << ERR_MSG << "RF Times must be in ascending order, starts with 0 and must not have duplicates values" << std::endl;
            return false;
        }
        // RF flip angles
        for(i=0; i<param.n_RF && ini.get("SCAN_PARAMETERS").has("RF_FA[" + std::to_string(i) + "]"); i++)
            param.RF_FA[i] = std::stof(ini.get("SCAN_PARAMETERS").get("RF_FA[" + std::to_string(i) + "]")) ;
        
        if(i != param.n_RF)
        {
            std::cout << ERR_MSG << "RF_FA and RF_ST must have the same number of elements" << std::endl;
            return false;
        }
        // RF phases
        for(i=0; i<param.n_RF && ini.get("SCAN_PARAMETERS").has("RF_PH[" + std::to_string(i) + "]"); i++)
            param.RF_PH[i] = std::stof(ini.get("SCAN_PARAMETERS").get("RF_PH[" + std::to_string(i) + "]")) ;

        if(i != param.n_RF)
        {
            std::cout << ERR_MSG << "RF_PH and RF_ST must have the same number of elements" << std::endl;
            return false;
        }
    }

    // ============== reading section SCAN_PARAMETERS ==============
    if(ini.has("STEADY_STATE"))
    {
        if(ini.get("STEADY_STATE").has("DUMMY_SCAN"))
            param.n_dummy_scan  = std::stoi(ini.get("STEADY_STATE").get("DUMMY_SCAN"));
        if(ini.get("STEADY_STATE").has("APPLY_FA/2"))
            param.enApplyFA2 = ini.get("STEADY_STATE").get("APPLY_FA/2").compare("0") != 0;
        if(ini.get("STEADY_STATE").has("PHASE_CYCLING"))
            param.phase_cycling = std::stof(ini.get("STEADY_STATE").get("PHASE_CYCLING")) ; // convert to radian
    }

    // ============== reading section SIMULATION_PARAMETERS ==============
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

    // ============== reading section TISSUE_PARAMETERS ==============
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


bool save_output(std::vector<float> &data, std::string output_filename, output_header hdr, std::vector<float> &additional_hdr)
{
    std::cout << "Saving output to: " << std::filesystem::absolute(output_filename) << std::endl;
    std::filesystem::path parent_path = std::filesystem::absolute(output_filename).parent_path();
    if (std::filesystem::is_directory(parent_path) == false)
    {
        std::cout << ERR_MSG << "cannot find directory " << parent_path.string() << ". Trying to create it." << std::endl;
        if(std::filesystem::create_directories(parent_path) == false)
        {
            std::cout << ERR_MSG << "cannot create directory " << parent_path.string() << std::endl;
            return false;
        }
    }

    std::ofstream file(output_filename, std::ios::out | std::ios::binary);
    if (file.is_open() == false)
    {
        std::cout << ERR_MSG << "cannot open file " << std::filesystem::absolute(output_filename) << std::endl;
        return false;
    }
    int32_t header_size = sizeof(output_header) + additional_hdr.size() * sizeof(additional_hdr[0]);
    file.write((char*)&header_size, sizeof(int32_t));
    file.write((char*)&hdr, sizeof(output_header));
    file.write((char*)additional_hdr.data(), additional_hdr.size() * sizeof(additional_hdr[0]));
    file.write((char*)data.data(), data.size() * sizeof(data[0]));
    file.close();
    return true; 
}

}
#endif  // __FILE_UTILS_H__