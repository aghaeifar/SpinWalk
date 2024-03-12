
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: file_utils.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#include <sstream>
#include <iterator>
#include <algorithm> 
#include <filesystem>
#include <boost/log/trivial.hpp>
#include "ini.h"
#include "file_utils.h"


//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
bool file_utils::read_config(std::string config_filename, simulation_parameters& param, std::vector<float>& sample_length_scales, std::map<std::string, std::vector<std::string> >& filenames)
{
    std::stringstream ss;
    if (std::filesystem::exists(config_filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Config-file does not exist: " << config_filename;
        return false;
    }
    // first, create a file instance
    mINI::INIFile file(config_filename);
    // next, create a structure that will hold data
    mINI::INIStructure ini;
    if (file.read(ini) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Problem reading config file: " << config_filename;
        return false;
    }

    BOOST_LOG_TRIVIAL(info) << "Reading config: " << config_filename;

    if(ini.has("parent"))
    {
        if(ini.get("parent").has("parent_config"))
        {
            std::filesystem::path parent_config(ini.get("parent").get("parent_config"));   
            if (parent_config.is_relative())
            {   // if parent_config is relative, make it absolute
                std::filesystem::path parent_path = std::filesystem::absolute(config_filename).parent_path();
                parent_config = parent_path / parent_config;
            }
            if (read_config(std::string(parent_config), param, sample_length_scales, filenames) == false)
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
            {               
                filenames.at("fieldmap").push_back(ini.get("files").get("fieldmap[" + std::to_string(i) + "]"));
                std::filesystem::path fieldmap(filenames.at("fieldmap").back());   
                if (fieldmap.is_relative())
                {   // if it is relative, make it absolute
                    std::filesystem::path parent_path = std::filesystem::absolute(config_filename).parent_path();
                    filenames.at("fieldmap").back() = (parent_path / fieldmap).string();
                }
            }
        }
        param.n_fieldmaps = filenames.at("fieldmap").size();

        // read all other files 
        for (std::map<std::string, std::vector<std::string>>::iterator it=filenames.begin(); it!=filenames.end(); ++it)
        {
            if( it->first.compare("fieldmap") == 0)
                continue; // fieldmap is already read
            it->second.assign(param.n_fieldmaps, ""); // this clear data from parent config files, is right?
            for (uint16_t i = 0; i<param.n_fieldmaps && ini.get("files").has(it->first + "[" + std::to_string(i) + "]"); i++)
            {
                it->second[i] = ini.get("files").get(it->first + "[" + std::to_string(i) + "]");
                if (std::all_of(it->second[i].begin(), it->second[i].end(), isspace))
                    break; // if it is empty, break the loop
                    
                std::filesystem::path f(it->second[i]);   
                if (f.is_relative())
                {   // if it is relative, make it absolute
                    std::filesystem::path parent_path = std::filesystem::absolute(config_filename).parent_path();
                    it->second[i] = (parent_path / f).string();
                }
            }
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
            ss.str(""); std::copy(param.TE, param.TE + param.n_TE, std::ostream_iterator<int>(ss, " "));
            BOOST_LOG_TRIVIAL(error) << "TE must be in ascending order and must not have duplicates or negative values: " << ss.str();
            return false;
        }        

        // ---------------- RF pulses (start times, Flip angles, phases and ) ----------------
        // RF start times
        for(i=0; i<MAX_RF && ini.get("SCAN_PARAMETERS").has("RF_ST[" + std::to_string(i) + "]"); i++)           
            param.RF_ST[i] = std::stof(ini.get("SCAN_PARAMETERS").get("RF_ST[" + std::to_string(i) + "]")) / param.dt;
        
        param.n_RF = i;
        // RF flip angles
        for(i=0; i<param.n_RF && ini.get("SCAN_PARAMETERS").has("RF_FA[" + std::to_string(i) + "]"); i++)
            param.RF_FA[i] = std::stof(ini.get("SCAN_PARAMETERS").get("RF_FA[" + std::to_string(i) + "]")) ;
        
        if(i != param.n_RF)
        {
            BOOST_LOG_TRIVIAL(error) << "RF_FA and RF_ST must have the same number of elements";
            return false;
        }
        // RF phases
        for(i=0; i<param.n_RF && ini.get("SCAN_PARAMETERS").has("RF_PH[" + std::to_string(i) + "]"); i++)
            param.RF_PH[i] = std::stof(ini.get("SCAN_PARAMETERS").get("RF_PH[" + std::to_string(i) + "]")) ;

        if(i != param.n_RF)
        {
            BOOST_LOG_TRIVIAL(error) << "RF_PH and RF_ST must have the same number of elements";
            return false;
        }

        // check RF start time conditions
        if (std::is_sorted(param.RF_ST, param.RF_ST + i) == false || 
            std::adjacent_find(param.RF_ST, param.RF_ST + i) != param.RF_ST + i || 
            param.RF_ST[0] != 0 || 
            param.n_RF == 0)
        {
            ss.str(""); std::copy(param.RF_ST, param.RF_ST + param.n_RF, std::ostream_iterator<int>(ss, " "));
            BOOST_LOG_TRIVIAL(error) << "RF Times must be in ascending order, starts with 0 and must not have duplicates values: " << ss.str();
            return false;
        }
        
        // ---------------- dephasing (start times, Flip angles ) ----------------
        // Dephase start times
        for(i=0; i<MAX_RF && ini.get("SCAN_PARAMETERS").has("DEPHASING_T[" + std::to_string(i) + "]"); i++)  
        {   // here we use try-catch since this field can be with empty value in the config file
            try{param.dephasing_T[i] = std::stof(ini.get("SCAN_PARAMETERS").get("DEPHASING_T[" + std::to_string(i) + "]")) / param.dt;} 
            catch(const std::exception& e){break;}
        }
 
        param.n_dephasing = i;   
        // Dephase flip angles
        for(i=0; i<param.n_dephasing && ini.get("SCAN_PARAMETERS").has("DEPHASING[" + std::to_string(i) + "]"); i++)
        {   
            try{param.dephasing[i] = std::stof(ini.get("SCAN_PARAMETERS").get("DEPHASING[" + std::to_string(i) + "]")) ;}
            catch(const std::exception& e){break;}
        }
        
        if(i != param.n_dephasing)
        {
            BOOST_LOG_TRIVIAL(error) << "DEPHASING and DEPHASING_T must have the same number of elements";
            return false;
        } 

        // check Dephase start time conditions
        if (std::is_sorted(param.dephasing_T, param.dephasing_T + i) == false || 
            std::adjacent_find(param.dephasing_T, param.dephasing_T + i) != param.dephasing_T + i )
        {
            ss.str(""); std::copy(param.dephasing_T, param.dephasing_T + param.n_dephasing, std::ostream_iterator<int>(ss, " ")); 
            BOOST_LOG_TRIVIAL(error) << "dephasing Times must be in ascending order and must not have duplicates values: " << ss.str();
            return false;
        } 
  
        // ---------------- Gradients (start times, strength (T/m) ) ----------------
        // Gradient start times
        for(i=0; i<MAX_GRADIENT && ini.get("SCAN_PARAMETERS").has("GRADIENT_T[" + std::to_string(i) + "]"); i++)  
        {   // here we use try-catch since this field can be with empty value in the config file   
            try{param.gradient_T[i] = std::stof(ini.get("SCAN_PARAMETERS").get("GRADIENT_T[" + std::to_string(i) + "]")) / param.dt; }
            catch(const std::exception& e){break;}  
        }

        param.n_gradient = i;
        // Gradient strength
        for(i=0; i<param.n_gradient && ini.get("SCAN_PARAMETERS").has("GRADIENT_XYZ[" + std::to_string(i) + "]"); i++)
        {
            try{
                std::istringstream iss(ini.get("SCAN_PARAMETERS").get("GRADIENT_XYZ[" + std::to_string(i) + "]"));
                int ind = i*3;
                while (iss >> *(param.gradient_xyz+ind++) && ind < (i+1)*3);
            }
            catch(const std::exception& e){break;}
        }
               
        if(i != param.n_gradient)
        {
            BOOST_LOG_TRIVIAL(error) << "GRADIENT_XYZ and GRADIENT_T must have the same number of elements";
            return false;
        } 

        // check Gradient start time conditions
        if (std::is_sorted(param.gradient_T, param.gradient_T + i) == false || 
            std::adjacent_find(param.gradient_T, param.gradient_T + i) != param.gradient_T + i )
        {
            ss.str(""); std::copy(param.gradient_T, param.gradient_T + param.n_gradient, std::ostream_iterator<int>(ss, " "));
            BOOST_LOG_TRIVIAL(error) << "Gradient Times must be in ascending order and must not have duplicates values: " << ss.str();
            return false;
        } 
        
    }

    // ============== reading section SCAN_PARAMETERS ==============
    if(ini.has("STEADY_STATE"))
    {
        if(ini.get("STEADY_STATE").has("DUMMY_SCAN"))
            param.n_dummy_scan  = std::stoi(ini.get("STEADY_STATE").get("DUMMY_SCAN"));
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
        if(ini.get("SIMULATION_PARAMETERS").has("CROSS_BOUNDARY"))
            param.enCrossBoundry = ini.get("SIMULATION_PARAMETERS").get("CROSS_BOUNDARY").compare("0") != 0;
        if(ini.get("SIMULATION_PARAMETERS").has("RECORD_TRAJECTORY"))
            param.enRecordTrajectory = ini.get("SIMULATION_PARAMETERS").get("RECORD_TRAJECTORY").compare("0") != 0;
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
        uint16_t i = 0;
        bool bT1exist = ini.get("TISSUE_PARAMETERS").has("T1[0]");
        param.n_T12 = bT1exist ? 0:param.n_T12;

        for(i=0; i<MAX_T12 && ini.get("TISSUE_PARAMETERS").has("T1[" + std::to_string(i) + "]"); i++)        
            param.T1[i] = std::stof(ini.get("TISSUE_PARAMETERS").get("T1[" + std::to_string(i) + "]"));

        param.n_T12 = bT1exist ? i:param.n_T12;

        for(i=0; i<MAX_T12 && ini.get("TISSUE_PARAMETERS").has("T2[" + std::to_string(i) + "]"); i++)        
            param.T2[i] = std::stof(ini.get("TISSUE_PARAMETERS").get("T2[" + std::to_string(i) + "]"));
        if(i != param.n_T12 && bT1exist)
        {
            BOOST_LOG_TRIVIAL(error) << "T1 and T2 must have the same number of elements (" << param.n_T12 << " vs " << i << ")";
            return false;
        }

        if(ini.get("TISSUE_PARAMETERS").has("MultiTissue"))
            param.enMultiTissue = ini.get("TISSUE_PARAMETERS").get("MultiTissue").compare("0") != 0;
    }   

    // ============== prep ==============
    param.prepare(); 

    if (param.n_dummy_scan < 0)
        param.n_dummy_scan = 5.0 * param.T1[0] / param.TR;

    return true;
}


bool file_utils::read_header(std::string filename, input_header &hdr_in)
{
    if(std::filesystem::exists(filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "File does not exist: " << filename;
        return false;
    }

    std::ifstream in_field(filename, std::ios::in | std::ios::binary);
    if (!in_field.is_open()) 
    {
        BOOST_LOG_TRIVIAL(error) << "Error reading header of " << filename;
        return false;
    }

    in_field.read((char*)&hdr_in, sizeof(input_header));
    return true;
}


bool file_utils::read_fieldmap(std::string fieldmap_filename, std::vector<float> &fieldmap, std::vector<uint8_t> &mask, simulation_parameters& param)
{
    if(std::filesystem::exists(fieldmap_filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Fieldmap file does not exist: " << fieldmap_filename;
        return false;
    }

    input_header hdr_in;
    BOOST_LOG_TRIVIAL(info) << "Opening fieldmap " << fieldmap_filename;
    std::ifstream in_field(fieldmap_filename, std::ios::in | std::ios::binary);
    if (!in_field.is_open()) 
    {
        BOOST_LOG_TRIVIAL(error) << "problem opening file " << fieldmap_filename;
        return false;
    }

    in_field.read((char*)&hdr_in, sizeof(input_header));

    BOOST_LOG_TRIVIAL(info) << "Size = " << hdr_in.fieldmap_size[0] << " x " << hdr_in.fieldmap_size[1] << " x " << hdr_in.fieldmap_size[2] << std::endl;
    BOOST_LOG_TRIVIAL(info) << "Length = " << hdr_in.sample_length[0]*1e6 << " x " << hdr_in.sample_length[1]*1e6 << " x " << hdr_in.sample_length[2]*1e6 << " um^3" << std::endl;
    
    std::copy(hdr_in.fieldmap_size, hdr_in.fieldmap_size + 3, param.fieldmap_size);
    std::copy(hdr_in.sample_length, hdr_in.sample_length + 3, param.sample_length);
    param.matrix_length = param.fieldmap_size[0] * param.fieldmap_size[1] * param.fieldmap_size[2];
    if (fieldmap.size() != param.matrix_length)
    {
        BOOST_LOG_TRIVIAL(info) << "Fieldmap size changed. Re-allocating memory...";
        BOOST_LOG_TRIVIAL(info) << "Old size: " << fieldmap.size();
        BOOST_LOG_TRIVIAL(info) << "New size: " << param.matrix_length;
        BOOST_LOG_TRIVIAL(info) << "New length (um): " << param.sample_length[0] * 1e6 << " " << param.sample_length[1] * 1e6 << " " << param.sample_length[2] * 1e6;
        fieldmap.resize(param.matrix_length);
        mask.resize(param.matrix_length);
    }

    BOOST_LOG_TRIVIAL(info) << "Reading...fieldmap...";
    in_field.read((char*)fieldmap.data(), sizeof(float) * param.matrix_length); BOOST_LOG_TRIVIAL(info) << "done...mask...";
    in_field.read((char*)mask.data(), sizeof(uint8_t) * param.matrix_length); BOOST_LOG_TRIVIAL(info) << "done.";
    in_field.close();
    return true;
}


bool file_utils::read_file(std::string filename, std::vector<float> &storage)
{
    if(std::filesystem::exists(filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "file does not exist: " << filename;
        return false;
    }

    BOOST_LOG_TRIVIAL(info) << "Reading " << filename;
    std::ifstream in_field(filename, std::ios::in | std::ios::binary);
    if (!in_field.is_open()) 
    {
        BOOST_LOG_TRIVIAL(error) << "error opening file " << filename;
        return false;
    }
    
    in_field.read((char*)storage.data(), sizeof(storage[0]) * storage.size());
    in_field.close();

    return true;
}


bool file_utils::save_output(std::vector<float> &data, std::string output_filename, output_header hdr, std::vector<float> &additional_hdr)
{
    BOOST_LOG_TRIVIAL(info) << "Saving output to: " << std::filesystem::absolute(output_filename);
    std::filesystem::path parent_path = std::filesystem::absolute(output_filename).parent_path();
    if (std::filesystem::is_directory(parent_path) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "cannot find directory " << parent_path.string() << ". Trying to create it.";
        if(std::filesystem::create_directories(parent_path) == false)
        {
            BOOST_LOG_TRIVIAL(error) << "cannot create directory " << parent_path.string();
            return false;
        }
    }

    std::ofstream file(output_filename, std::ios::out | std::ios::binary);
    if (file.is_open() == false)
    {
        BOOST_LOG_TRIVIAL(error) << "cannot open file " << std::filesystem::absolute(output_filename);
        return false;
    }

    int32_t add_hdr_size = additional_hdr.size() * sizeof(additional_hdr[0]);
    int32_t header_size  = sizeof(output_header) + add_hdr_size;
    file.write((char*)&header_size, sizeof(int32_t));
    file.write((char*)&hdr, sizeof(output_header));
    file.write((char*)additional_hdr.data(), add_hdr_size);
    file.write((char*)data.data(), data.size() * sizeof(data[0]));
    file.close();
    return true; 
}
