
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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "file_utils.h"

uint8_t find_max(const std::vector<uint8_t> &data);
//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
std::filesystem::path output_dir("./output");

bool file_utils::read_config(std::string config_filename, simulation_parameters& param, std::vector<double>& sample_length_scales, std::map<std::string, std::vector<std::string> >& filenames)
{
    std::stringstream ss;
    if (std::filesystem::exists(config_filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Config-file does not exist: " << config_filename;
        return false;
    }
    
    std::string cf_name = std::filesystem::path(config_filename).filename();
    BOOST_LOG_TRIVIAL(info) << "Reading config: " << config_filename;
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(config_filename, pt);

    if(pt.get_child_optional("PARENT.PARENT_CONFIG"))
    {
        std::filesystem::path parent_config(pt.get<std::string>("PARENT.PARENT_CONFIG", ""));   
        if (parent_config.empty() == false)
        {
            if (parent_config.is_relative())
            {   // if parent_config is relative, make it absolute
                std::filesystem::path parent_path = std::filesystem::absolute(config_filename).parent_path();
                parent_config = parent_path / parent_config;
            }
            if (read_config(parent_config.string(), param, sample_length_scales, filenames) == false)
                return false;
        }
    }

    std::string seq_name = pt.get<std::string>("SEQ_NAME", "");
    // ============== reading section FILES ==============
    std::vector<std::string> file_paths;
    for (const auto& str : {std::string("FIELDMAP"), std::string("XYZ0"), std::string("M0")}) 
    {
        file_paths.clear();
        for (uint16_t i = 0; pt.get_child_optional("FILES." + str + "[" + std::to_string(i) + "]") ; i++) 
            file_paths.push_back(pt.get<std::string>("FILES." + str + "[" + std::to_string(i) + "]", ""));
        // remove empty strings
        file_paths.erase(std::remove(file_paths.begin(), file_paths.end(), ""), file_paths.end());
        // if it is relative, make it absolute
        for (auto& path : file_paths) 
            if (std::filesystem::path(path).is_relative()) 
                path = (std::filesystem::absolute(config_filename).parent_path() / path).string();
        // replace if there is any
        if (file_paths.size() > 0)
            filenames[str] = file_paths;
    }
    param.n_fieldmaps = filenames["FIELDMAP"].size();
      
    // output directory 
    output_dir = std::filesystem::path(pt.get<std::string>("FILES.OUTPUT_DIR", output_dir.string()));
    if (output_dir.is_relative())
        output_dir = std::filesystem::absolute(config_filename).parent_path() / output_dir;
    // generae names for m1
    file_paths.clear();
    for (const auto& path : filenames["FIELDMAP"])
        file_paths.push_back((output_dir / (seq_name + "_m1_" + std::filesystem::path(path).filename().string())).string());
    filenames["M1"] = file_paths;
    // generae names for xyz1  
    file_paths.clear();
    for (const auto& path : filenames["FIELDMAP"])
        file_paths.push_back((output_dir / (seq_name + "_xyz1_" + std::filesystem::path(path).filename().string())).string());
    filenames["XYZ1"] = file_paths;
    // generae names for T  
    file_paths.clear();
    for (const auto& path : filenames["FIELDMAP"])
        file_paths.push_back((output_dir / (seq_name + "_T_" + std::filesystem::path(path).filename().string())).string());
    filenames["T"] = file_paths;
    
    // check header of fieldmap matches
    file_utils::input_header hdr_in;
    if(file_utils::read_header(filenames.at("FIELDMAP")[0], hdr_in) == false)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "reading header of fieldmap " << filenames.at("FIELDMAP")[0] << " failed. Aborting...!";
        return false;
    }
    std::copy(hdr_in.fieldmap_size, hdr_in.fieldmap_size+3, param.fieldmap_size);
    std::copy(hdr_in.sample_length, hdr_in.sample_length+3, param.sample_length);
    param.file_size = hdr_in.file_size;
    param.fieldmap_exist = hdr_in.fieldmap_exists;
    param.mask_exist = hdr_in.mask_exists;

    // ============== reading section SCAN_PARAMETERS ==============
    param.TR = pt.get("SCAN_PARAMETERS.TR", param.TR);
    param.dt = pt.get("SCAN_PARAMETERS.DWELL_TIME", param.dt);

    uint16_t i=0, j=0;
    // ---------------- Echo times ----------------       
    for(i=0; i<MAX_TE && pt.get("SCAN_PARAMETERS.TE[" + std::to_string(i) + "]", -1.f) != -1.f; i++)        
        param.TE[i] = pt.get<double>("SCAN_PARAMETERS.TE[" + std::to_string(i) + "]", 0.) / param.dt;
    param.n_TE = i==0?param.n_TE:i;
    // check TE conditions
    if (std::is_sorted(param.TE, param.TE + i) == false || std::adjacent_find(param.TE, param.TE + i) != param.TE + i || param.TE[0] < 0 || param.n_TE == 0)
    {
        ss.str(""); std::copy(param.TE, param.TE + param.n_TE, std::ostream_iterator<int>(ss, " "));
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "TE must be in ascending order and must not have duplicates or negative values: " << ss.str();
        return false;
    }     

    // ---------------- RF pulses (start times, Flip angles, phases and ) ----------------
    // RF start times
    for(i=0; i<MAX_RF && pt.get("SCAN_PARAMETERS.RF_ST[" + std::to_string(i) + "]", -1.f) != -1.f; i++)           
        param.RF_ST[i] = pt.get<double>("SCAN_PARAMETERS.RF_ST[" + std::to_string(i) + "]", 0.) / param.dt;
    param.n_RF = i==0?param.n_RF:i;
    // RF flip angles    
    for(j=0; j<param.n_RF && pt.get("SCAN_PARAMETERS.RF_FA[" + std::to_string(j) + "]", -1.f) != -1.f; j++)
        param.RF_FA[j] = pt.get<double>("SCAN_PARAMETERS.RF_FA[" + std::to_string(j) + "]", 0.) ; 
    if(j !=i && j != param.n_RF)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "RF_FA and RF_ST must have the same number of elements " << j << " vs " << param.n_RF;
        return false;
    }
    // RF phases
    for(j=0; j<param.n_RF && pt.get("SCAN_PARAMETERS.RF_PH[" + std::to_string(j) + "]", -1.f) != -1.f; j++)
        param.RF_PH[j] = pt.get<double>("SCAN_PARAMETERS.RF_PH[" + std::to_string(j) + "]", 0.) ;
    if(j !=i && j != param.n_RF)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "RF_PH and RF_ST must have the same number of elements " << j << " vs " << param.n_RF; 
        return false;
    }
    // check RF start time conditions
    if (std::is_sorted(param.RF_ST, param.RF_ST + i) == false || std::adjacent_find(param.RF_ST, param.RF_ST + i) != param.RF_ST + i || param.RF_ST[0] != 0 || param.n_RF == 0)
    {
        ss.str(""); std::copy(param.RF_ST, param.RF_ST + param.n_RF, std::ostream_iterator<int>(ss, " "));
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "RF Times must be in ascending order, starts with 0 and must not have duplicates values: " << ss.str();
        return false;
    }
 
    // ---------------- dephasing (start times, Flip angles ) ----------------
    // Dephase start times
    for(i=0; i<MAX_RF && pt.get("SCAN_PARAMETERS.DEPHASING_T[" + std::to_string(i) + "]", -1.f) != -1.f; i++)  
        param.dephasing_T[i] = pt.get<double>("SCAN_PARAMETERS.DEPHASING_T[" + std::to_string(i) + "]", 0.) / param.dt;
    param.n_dephasing = i==0?param.n_dephasing:i;
    // Dephase flip angles
    for(j=0; j<param.n_dephasing && pt.get_child_optional("SCAN_PARAMETERS.DEPHASING[" + std::to_string(j) + "]"); j++)
        param.dephasing[j] = pt.get<double>("SCAN_PARAMETERS.DEPHASING[" + std::to_string(j) + "]", 0.) ;
    if(j !=i && j != param.n_dephasing)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "DEPHASING and DEPHASING_T must have the same number of elements " << j << " vs " << param.n_dephasing;
        return false;
    }
    // check Dephase start time conditions
    if (std::is_sorted(param.dephasing_T, param.dephasing_T + i) == false || std::adjacent_find(param.dephasing_T, param.dephasing_T + i) != param.dephasing_T + i)
    {
        ss.str(""); std::copy(param.dephasing_T, param.dephasing_T + param.n_dephasing, std::ostream_iterator<int>(ss, " ")); 
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "dephasing Times must be in ascending order and must not have duplicates values: " << ss.str();
        return false;
    }

    // ---------------- Gradients (start times, strength (T/m) ) ----------------
    // Gradient start times
    for(i=0; i<MAX_GRADIENT && pt.get("SCAN_PARAMETERS.GRADIENT_T[" + std::to_string(i) + "]", -1.f) != -1.f; i++)  
        param.gradient_T[i] = pt.get<double>("SCAN_PARAMETERS.GRADIENT_T[" + std::to_string(i) + "]", 0.) / param.dt;
    param.n_gradient = i==0?param.n_gradient:i;
    // Gradient strength
    for(j=0; j<param.n_gradient && pt.get_child_optional("SCAN_PARAMETERS.GRADIENT_XYZ[" + std::to_string(j) + "]"); j++)
    {
        std::istringstream iss(pt.get<std::string>("SCAN_PARAMETERS.GRADIENT_XYZ[" + std::to_string(j) + "]", "0.0 0.0 0.0"));
        int ind = j*3;
        while (iss >> *(param.gradient_xyz+ind++) && ind < (j+1)*3);
    }
    if(j !=i && j != param.n_gradient)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "GRADIENT_XYZ and GRADIENT_T must have the same number of elements " << j << " vs " << param.n_gradient;
        return false;
    }
    // check Gradient start time conditions
    if (std::is_sorted(param.gradient_T, param.gradient_T + i) == false || std::adjacent_find(param.gradient_T, param.gradient_T + i) != param.gradient_T + i )
    {
        ss.str(""); std::copy(param.gradient_T, param.gradient_T + param.n_gradient, std::ostream_iterator<int>(ss, " ")); 
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "Gradient Times must be in ascending order and must not have duplicates values: " << ss.str();
        return false;
    }

    // ============== reading section SCAN_PARAMETERS ==============
    param.n_dummy_scan = pt.get<int32_t>("STEADY_STATE.DUMMY_SCAN", param.n_dummy_scan);
    param.phase_cycling = pt.get("STEADY_STATE.PHASE_CYCLING", param.phase_cycling); 

    // ============== reading section SIMULATION_PARAMETERS ==============
    param.B0                    = pt.get("SIMULATION_PARAMETERS.B0", param.B0);
    param.seed                  = pt.get<float>("SIMULATION_PARAMETERS.SEED", param.seed);
    param.n_spins               = pt.get<float>("SIMULATION_PARAMETERS.NUMBER_OF_SPINS", param.n_spins); // template type must be float since input can be of form scientific notation
    param.enCrossFOV            = pt.get("SIMULATION_PARAMETERS.CROSS_FOV", param.enCrossFOV);
    param.enRecordTrajectory    = pt.get("SIMULATION_PARAMETERS.RECORD_TRAJECTORY", param.enRecordTrajectory);
    param.diffusion_const       = pt.get("SIMULATION_PARAMETERS.DIFFUSION_CONSTANT", param.diffusion_const);
    param.max_iterations        = pt.get<float>("SIMULATION_PARAMETERS.MAX_ITERATIONS", param.max_iterations);

    std::vector<double> sls;
    for(i=0; i<MAX_RF && pt.get_child_optional("SIMULATION_PARAMETERS.SAMPLE_LENGTH_SCALES[" + std::to_string(i) + "]"); i++) 
        sls.push_back(pt.get("SIMULATION_PARAMETERS.SAMPLE_LENGTH_SCALES[" + std::to_string(i) + "]", 0.));

    if(sls.size() > 0)
        sample_length_scales = sls;
    param.n_sample_length_scales = sample_length_scales.size();

    // ============== reading section TISSUE_PARAMETERS ==============
    // T1 & T2
    for(i=0; i<MAX_T12 && pt.get_child_optional("TISSUE_PARAMETERS.T1[" + std::to_string(i) + "]"); i++)  
        param.T1[i] = pt.get("TISSUE_PARAMETERS.T1[" + std::to_string(i) + "]", 10000.f);
    param.n_T12 = i==0?param.n_T12:i;
    for(j=0; j<i && pt.get_child_optional("TISSUE_PARAMETERS.T2[" + std::to_string(j) + "]"); j++)  
        param.T2[j] = pt.get("TISSUE_PARAMETERS.T2[" + std::to_string(j) + "]", 10000.f);

    if(j != i && j!=param.n_T12)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "T1 and T2 must have the same number of elements (" << j << " vs " << param.n_T12 << ")";
        return false;
    }

    // Cross Tissue Probability
    ss.str("");
    for(i=0; i<MAX_TISSUE_TYPE && pt.get<std::string>("TISSUE_PARAMETERS.P_XY[" + std::to_string(i) + "]", std::string("?")) != "?"; i++)
        ss << pt.get<std::string>("TISSUE_PARAMETERS.P_XY[" + std::to_string(i) + "]") << " ";
    param.n_tissue_type = i==0?param.n_tissue_type:i;
    for(i=0; i<param.n_tissue_type; i++)
        for(j=0; j<param.n_tissue_type; j++)
            ss >> param.pXY[j+i*param.n_tissue_type];
    
    // ============== prep ==============
    param.prepare(); 

    return true;
}


bool file_utils::read_header(std::string filename, input_header &hdr_in)
{
    if(std::filesystem::exists(filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "File does not exist: " << filename;
        return false;
    }

    hdr_in.file_size = std::filesystem::file_size(filename);

    std::ifstream in_field(filename, std::ios::in | std::ios::binary);
    if (!in_field.is_open()) 
    {
        BOOST_LOG_TRIVIAL(error) << "Error reading header of " << filename;
        return false;
    }
    float buff;
    in_field.read((char*)hdr_in.fieldmap_size, sizeof(hdr_in.fieldmap_size));
    in_field.read((char*)&buff, sizeof(buff)); hdr_in.sample_length[0] = buff;
    in_field.read((char*)&buff, sizeof(buff)); hdr_in.sample_length[1] = buff;
    in_field.read((char*)&buff, sizeof(buff)); hdr_in.sample_length[2] = buff;

    size_t matrix_length = hdr_in.fieldmap_size[0] * hdr_in.fieldmap_size[1] * hdr_in.fieldmap_size[2];
    size_t header_size = 6 * sizeof(float); // this should be later corrected to match input_header
    if (hdr_in.file_size - header_size == matrix_length * sizeof(float))
    {
        BOOST_LOG_TRIVIAL(error) << "Mask with labeled tissues is not provided in " << filename;
        hdr_in.mask_exists = false;
    }
    if (hdr_in.file_size - header_size == matrix_length)
    {
        BOOST_LOG_TRIVIAL(warning) << "Off-resonance map is not provided in " << filename;
        BOOST_LOG_TRIVIAL(warning) << "Simulation will be performed with a perfect homogeneous field.";
        hdr_in.fieldmap_exists = false;
    }

    return true;
}


bool file_utils::read_fieldmap(std::string fieldmap_filename, std::vector<float> &fieldmap, std::vector<uint8_t> &mask, simulation_parameters& param)
{
    if(std::filesystem::exists(fieldmap_filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Fieldmap file does not exist: " << fieldmap_filename;
        return false;
    }

    size_t file_size = std::filesystem::file_size(fieldmap_filename);
    if (file_size != param.file_size)
    {
        BOOST_LOG_TRIVIAL(error) << "All field map files must be of same size: " << file_size << " vs " << param.file_size;
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

    float buff;
    in_field.read((char*)hdr_in.fieldmap_size, sizeof(hdr_in.fieldmap_size));
    in_field.read((char*)&buff, sizeof(buff)); hdr_in.sample_length[0] = buff;
    in_field.read((char*)&buff, sizeof(buff)); hdr_in.sample_length[1] = buff;
    in_field.read((char*)&buff, sizeof(buff)); hdr_in.sample_length[2] = buff;

    BOOST_LOG_TRIVIAL(info) << "Size = " << hdr_in.fieldmap_size[0] << " x " << hdr_in.fieldmap_size[1] << " x " << hdr_in.fieldmap_size[2] << std::endl;
    BOOST_LOG_TRIVIAL(info) << "Length = " << hdr_in.sample_length[0]*1e6 << " x " << hdr_in.sample_length[1]*1e6 << " x " << hdr_in.sample_length[2]*1e6 << " um^3" << std::endl;

    if (param.fieldmap_exist)
    {
        BOOST_LOG_TRIVIAL(info) << "Reading...fieldmap...";
        in_field.read((char*)fieldmap.data(), sizeof(fieldmap[0]) * param.matrix_length); 
    }

    if (param.mask_exist == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Mask must exist in the fieldmap file: " << fieldmap_filename << ". Aborting...!";
        return false;
    }
    BOOST_LOG_TRIVIAL(info) << "Reading...mask...";
    in_field.read((char*)mask.data(), sizeof(mask[0]) * param.matrix_length);
    in_field.close();

    // int n_tissue = *std::max_element(mask.begin(), mask.end(),  [](const uint8_t &x,const uint8_t &y) {return x<y;}) + 1;
    int n_tissue = find_max(mask) + 1;
    if (n_tissue > param.n_tissue_type)
    {
        BOOST_LOG_TRIVIAL(error) << "The number of tissue types in the mask does not match the number of tissue types in the config file: " << n_tissue << " vs " << param.n_tissue_type;
        return false;
    }
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


bool file_utils::save_output(char *data, size_t bytes, std::string output_filename, output_header hdr, std::vector<double> &additional_hdr)
{
    BOOST_LOG_TRIVIAL(info) << "Saving output to: " << std::filesystem::absolute(output_filename);
    std::filesystem::path parent_path = std::filesystem::absolute(output_filename).parent_path();
    if (std::filesystem::is_directory(parent_path) == false)
    {
        BOOST_LOG_TRIVIAL(warning) << "cannot find directory " << parent_path.string() << ". Trying to create it.";
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
    file.write(data, bytes);
    file.close();
    return true; 
}
