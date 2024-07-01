
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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>
#include "file_utils.h"
#include <boost/log/trivial.hpp> 
#include <highfive/highfive.hpp>
#include "simulation_parameters.h"

uint8_t find_max(const std::vector<uint8_t> &data);
//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
std::filesystem::path output_dir("./output");

bool file_utils::read_config(std::string config_filename, simulation_parameters *param, std::vector<double>& fov_scale, std::map<std::string, std::vector<std::string> >& filenames, bool isParentConfig)
{
    std::stringstream ss;
    if (std::filesystem::exists(config_filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Config-file does not exist: " << config_filename;
        return false;
    }
    
    std::string cf_name = std::filesystem::path(config_filename).filename().string();
    BOOST_LOG_TRIVIAL(info) << "Reading config: " << config_filename;
    boost::property_tree::ptree pt;

    try{ 
        boost::property_tree::ini_parser::read_ini(config_filename, pt);
    }catch (std::exception& e) {
        std::cout << ERR_MSG << e.what() << "\n";
        return false;
    } 

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
            if (read_config(parent_config.string(), param, fov_scale, filenames, true) == false)
                return false;

        BOOST_LOG_TRIVIAL(info) << "Back to reading config: " << config_filename;
        }
    }

    std::string seq_name = pt.get<std::string>("SEQ_NAME", "");
    // ============== reading section FILES ==============
    std::vector<std::string> file_paths;
    for (const auto& str : {std::string("PHANTOM"), std::string("XYZ0"), std::string("M0")}) 
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
            filenames[boost::algorithm::to_lower_copy(str)] = file_paths;
    }
    param->n_fieldmaps = filenames["phantom"].size();
    if(isParentConfig == false && param->n_fieldmaps == 0)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "No fieldmap is provided. Aborting...!";
        return false;
    }
      
    // output directory 
    output_dir = std::filesystem::path(pt.get<std::string>("FILES.OUTPUT_DIR", output_dir.string()));
    if (output_dir.is_relative())
        output_dir = std::filesystem::absolute(config_filename).parent_path() / output_dir;
    // generate names for output
    file_paths.clear();
    for (const auto& path : filenames["phantom"])
    {
        auto f = output_dir / (seq_name + "_" + std::filesystem::path(path).filename().string());
        f.replace_extension(".h5");
        file_paths.push_back(f.string());
    }
    filenames["output"] = file_paths;
     
    if(param->n_fieldmaps > 0)
    {
        // check header of fieldmap matches
        auto dims = file_utils::get_size_h5(filenames.at("phantom")[0], "mask");
        if(dims.size() != 3)
        {
            BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "mask must be 3D but is " << dims.size() << "D filename:"<< filenames.at("phantom")[0] << ". Aborting...!";
            return false;
        }
        std::vector<float> fov(3, 0);
        file_utils::read_h5(filenames.at("phantom")[0], fov.data(), "fov", "float");
        std::copy(dims.begin(), dims.end(), param->fieldmap_size);
        std::copy(fov.begin(), fov.end(), param->fov);

        param->fieldmap_exist = file_utils::get_size_h5(filenames.at("phantom")[0], "fieldmap").size() == 3;
        param->mask_exist     = file_utils::get_size_h5(filenames.at("phantom")[0], "mask").size() == 3;
    }

    // ============== reading section SCAN_PARAMETERS ==============
    param->TR_us = (int32_t)pt.get<float>("SCAN_PARAMETERS.TR", param->TR_us);
    param->timestep_us = pt.get("SCAN_PARAMETERS.TIME_STEP", param->timestep_us);

    uint16_t i=0, j=0;
    // ---------------- Echo times ----------------       
    for(i=0; i<MAX_TE && pt.get<float>("SCAN_PARAMETERS.TE[" + std::to_string(i) + "]", -1) >= 0 ; i++)        
        param->TE_us[i] = (int32_t)pt.get<float>("SCAN_PARAMETERS.TE[" + std::to_string(i) + "]", 0) / param->timestep_us;
      
    param->n_TE = (i==0?param->n_TE:i);
    // check TE conditions
    if (std::is_sorted(param->TE_us, param->TE_us + param->n_TE) == false || 
        std::adjacent_find(param->TE_us, param->TE_us + param->n_TE) != param->TE_us + param->n_TE || 
        param->TE_us[0] < 0 || param->n_TE == 0)
    {
        ss.str(""); std::copy(param->TE_us, param->TE_us + param->n_TE, std::ostream_iterator<int>(ss, " "));
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "TE must be in ascending order and must not have duplicates or negative values: " << ss.str();
        return false;
    }     

    // ---------------- RF pulses (start times, Flip angles, phases and ) ----------------
    // RF start times
    for(i=0; i<MAX_RF && pt.get<float>("SCAN_PARAMETERS.RF_T[" + std::to_string(i) + "]", -1) >= 0; i++)           
        param->RF_us[i] = (int32_t)pt.get<float>("SCAN_PARAMETERS.RF_T[" + std::to_string(i) + "]", 0.) / param->timestep_us;
    param->n_RF = i==0?param->n_RF:i;
    // RF flip angles    
    for(j=0; j<param->n_RF && pt.get("SCAN_PARAMETERS.RF_FA[" + std::to_string(j) + "]", -1.f) != -1.f; j++)
        param->RF_FA_deg[j] = pt.get<double>("SCAN_PARAMETERS.RF_FA[" + std::to_string(j) + "]", 0.) ; 
    if(j !=i && j != param->n_RF)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "RF_FA and RF_us must have the same number of elements " << j << " vs " << param->n_RF;
        return false;
    }
    // RF phases
    for(j=0; j<param->n_RF && pt.get("SCAN_PARAMETERS.RF_PH[" + std::to_string(j) + "]", -1.f) != -1.f; j++)
        param->RF_PH_deg[j] = pt.get<double>("SCAN_PARAMETERS.RF_PH[" + std::to_string(j) + "]", 0.) ;
    if(j !=i && j != param->n_RF)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "RF_PH and RF_us must have the same number of elements " << j << " vs " << param->n_RF; 
        return false;
    }
    // check RF start time conditions
    if (std::is_sorted(param->RF_us, param->RF_us + i) == false || 
        std::adjacent_find(param->RF_us, param->RF_us + i) != param->RF_us + i || 
        param->RF_us[0] != 0 || param->n_RF == 0)
    {
        ss.str(""); std::copy(param->RF_us, param->RF_us + param->n_RF, std::ostream_iterator<int>(ss, " "));
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "RF Times must be in ascending order, starts with 0 and must not have duplicates values: " << ss.str();
        return false;
    }
 
    // ---------------- dephasing (start times, Flip angles ) ----------------
    // Dephase start times
    for(i=0; i<MAX_RF && pt.get<float>("SCAN_PARAMETERS.DEPHASING_T[" + std::to_string(i) + "]", -1.f) != -1.f; i++)  
        param->dephasing_us[i] = (int32_t)pt.get<float>("SCAN_PARAMETERS.DEPHASING_T[" + std::to_string(i) + "]", 0.) / param->timestep_us;
    param->n_dephasing = i==0?param->n_dephasing:i;
    // Dephase flip angles
    for(j=0; j<param->n_dephasing && pt.get_child_optional("SCAN_PARAMETERS.DEPHASING[" + std::to_string(j) + "]"); j++)
        param->dephasing_deg[j] = pt.get<double>("SCAN_PARAMETERS.DEPHASING[" + std::to_string(j) + "]", 0.) ;
    if(j !=i && j != param->n_dephasing)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "DEPHASING and DEPHASING_T must have the same number of elements " << j << " vs " << param->n_dephasing;
        return false;
    }
    // check Dephase start time conditions
    if (std::is_sorted(param->dephasing_us, param->dephasing_us + i) == false || 
        std::adjacent_find(param->dephasing_us, param->dephasing_us + i) != param->dephasing_us + i)
    {
        ss.str(""); std::copy(param->dephasing_us, param->dephasing_us + param->n_dephasing, std::ostream_iterator<int>(ss, " ")); 
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "dephasing Times must be in ascending order and must not have duplicates values: " << ss.str();
        return false;
    }

    // ---------------- Gradients (start times, strength (T/m) ) ----------------
    // Gradient start times
    for(i=0; i<MAX_GRADIENT && pt.get<float>("SCAN_PARAMETERS.GRADIENT_T[" + std::to_string(i) + "]", -1.f) != -1.f; i++)  
        param->gradient_us[i] = (int32_t)pt.get<float>("SCAN_PARAMETERS.GRADIENT_T[" + std::to_string(i) + "]", 0.) / param->timestep_us;
    param->n_gradient = i==0?param->n_gradient:i;
    // Gradient strength
    for(j=0; j<param->n_gradient && pt.get_child_optional("SCAN_PARAMETERS.GRADIENT_XYZ[" + std::to_string(j) + "]"); j++)
    {
        std::istringstream iss(pt.get<std::string>("SCAN_PARAMETERS.GRADIENT_XYZ[" + std::to_string(j) + "]", "0.0 0.0 0.0"));
        int ind = j*3;
        while (iss >> *(param->gradient_mTm+ind++) && ind < (j+1)*3);
    }
    if(j !=i && j != param->n_gradient)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "GRADIENT_XYZ and GRADIENT_T must have the same number of elements " << j << " vs " << param->n_gradient;
        return false;
    }
    // check Gradient start time conditions
    if (std::is_sorted(param->gradient_us, param->gradient_us + i) == false ||
        std::adjacent_find(param->gradient_us, param->gradient_us + i) != param->gradient_us + i )
    {
        ss.str(""); std::copy(param->gradient_us, param->gradient_us + param->n_gradient, std::ostream_iterator<int>(ss, " ")); 
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "Gradient Times must be in a strickly ascending order and must not have duplicates values: " << ss.str();
        return false;
    }

    // ============== reading section SCAN_PARAMETERS ==============
    param->n_dummy_scan = (int32_t)pt.get<float>("SCAN_PARAMETERS.DUMMY_SCAN", param->n_dummy_scan);
    param->phase_cycling = pt.get("SCAN_PARAMETERS.PHASE_CYCLING", param->phase_cycling); 

    // ============== reading section SIMULATION_PARAMETERS ==============
    param->B0                    = pt.get("SIMULATION_PARAMETERS.B0", param->B0);
    param->seed                  = pt.get<float>("SIMULATION_PARAMETERS.SEED", param->seed);
    param->n_spins               = pt.get<float>("SIMULATION_PARAMETERS.NUMBER_OF_SPINS", param->n_spins); // template type must be float since input can be of form scientific notation
    param->enCrossFOV            = pt.get("SIMULATION_PARAMETERS.CROSS_FOV", param->enCrossFOV);
    param->enRecordTrajectory    = pt.get("SIMULATION_PARAMETERS.RECORD_TRAJECTORY", param->enRecordTrajectory);
    param->enProfiling           = pt.get("SIMULATION_PARAMETERS.PROFILING", param->enProfiling);
    param->max_iterations        = pt.get<float>("SIMULATION_PARAMETERS.MAX_ITERATIONS", param->max_iterations);

    std::vector<double> sls;
    for(i=0; i<MAX_RF && pt.get_child_optional("SIMULATION_PARAMETERS.FOV_SCALE[" + std::to_string(i) + "]"); i++) 
        sls.push_back(pt.get("SIMULATION_PARAMETERS.FOV_SCALE[" + std::to_string(i) + "]", 0.));

    if(sls.size() > 0)
        fov_scale = sls;
    param->n_fov_scale = fov_scale.size();

    // ============== reading section TISSUE_PARAMETERS ==============
    // Diffusivity
    for(i=0; i<MAX_TISSUE_TYPE && pt.get_child_optional("TISSUE_PARAMETERS.DIFFUSIVITY[" + std::to_string(i) + "]"); i++)  
        param->diffusivity[i] = pt.get<double>("TISSUE_PARAMETERS.DIFFUSIVITY[" + std::to_string(i) + "]", param->diffusivity[i]);
    param->n_tissue_type = i==0?param->n_tissue_type:i;
    // T1 & T2
    for(i=0; i<param->n_tissue_type && pt.get_child_optional("TISSUE_PARAMETERS.T1[" + std::to_string(i) + "]"); i++)  
        param->T1_ms[i] = pt.get<float>("TISSUE_PARAMETERS.T1[" + std::to_string(i) + "]", 10000);
    if (i != 0 && i != param->n_tissue_type)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "T1 and diffusivity must have the same number of elements (" << i << " vs " << param->n_tissue_type << ")";
        return false;
    }

    for(i=0; i<param->n_tissue_type && pt.get_child_optional("TISSUE_PARAMETERS.T2[" + std::to_string(i) + "]"); i++)  
        param->T2_ms[i] = pt.get<float>("TISSUE_PARAMETERS.T2[" + std::to_string(i) + "]", 10000);
    if (i != 0 && i != param->n_tissue_type)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "T2 and diffusivity must have the same number of elements (" << i << " vs " << param->n_tissue_type << ")";
        return false;
    }

    // Cross Tissue Probability
    ss.str("");
    for(i=0; i<MAX_TISSUE_TYPE && pt.get<std::string>("TISSUE_PARAMETERS.P_XY[" + std::to_string(i) + "]", std::string("?")) != "?"; i++)
        ss << pt.get<std::string>("TISSUE_PARAMETERS.P_XY[" + std::to_string(i) + "]") << " ";
    if (i != 0 && i != param->n_tissue_type)
    {
        BOOST_LOG_TRIVIAL(error) << cf_name << ") " << "P_XY must have " << param->n_tissue_type << " rows (and columns) but has " << i << " rows";
        return false;
    }
    for(i=0; i<param->n_tissue_type; i++)
        for(j=0; j<param->n_tissue_type; j++)
            ss >> param->pXY[j+i*param->n_tissue_type];
    
    // ============== prep ==============
    if(isParentConfig == false)
        return param->prepare(); 
    else
        return true;
}


bool file_utils::read_phantom(std::string phantom_filename, std::vector<float> &fieldmap, std::vector<uint8_t> &mask, simulation_parameters *param)
{
    if(std::filesystem::exists(phantom_filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Fieldmap file does not exist: " << phantom_filename;
        return false;
    }

    BOOST_LOG_TRIVIAL(info) << "Opening fieldmap " << phantom_filename;
    std::vector<size_t> dims;
    if (param->fieldmap_exist)
    {
        BOOST_LOG_TRIVIAL(info) << "Reading...fieldmap...";
        dims = get_size_h5(phantom_filename, "fieldmap");
        if (product(dims) != fieldmap.size())
        {
            BOOST_LOG_TRIVIAL(warning) << "Phantom size has changed in " << phantom_filename;
            BOOST_LOG_TRIVIAL(warning) << "Fieldmap size does not match the expected size: " << product(dims) << " vs " << fieldmap.size() << ". resize it...!";
            fieldmap.resize(product(dims));
        }
        read_h5(phantom_filename, fieldmap.data(), "fieldmap", "float");
    }

    if (param->mask_exist)
    {
        BOOST_LOG_TRIVIAL(info) << "Reading...mask...";
        dims = get_size_h5(phantom_filename, "mask");
        if (product(dims) != mask.size())
        {
            BOOST_LOG_TRIVIAL(warning) << "Hint in reading mask: " << phantom_filename;
            BOOST_LOG_TRIVIAL(warning) << "Mask size does not match the expected size: " << product(dims) << " vs " << mask.size() << ". Aborting...!";
            mask.resize(product(dims));
        }
        read_h5(phantom_filename, mask.data(), "mask", "uint8_t");
    }

    if (param->mask_exist && param->fieldmap_exist && mask.size() != fieldmap.size())
    {
        BOOST_LOG_TRIVIAL(error) << "Fieldmap and mask sizes do not match: " << mask.size() << " vs " << fieldmap.size();
        return false;
    }    

    std::vector<float> fov(3, 0);
    file_utils::read_h5(phantom_filename, fov.data(), "fov", "float");
    std::copy(fov.begin() , fov.end() , param->fov);
    std::copy(dims.begin(), dims.end(), param->fieldmap_size);
    BOOST_LOG_TRIVIAL(info) << "Size = " << dims[0] << " x " << dims[1] << " x " << dims[2] << std::endl;
    BOOST_LOG_TRIVIAL(info) << "FoV = " << param->fov[0]*1e6 << " x " << param->fov[1]*1e6 << " x " << param->fov[2]*1e6 << " um^3" << std::endl;

    int n_tissue = find_max(mask) + 1;
    if (n_tissue > param->n_tissue_type)
    {
        BOOST_LOG_TRIVIAL(error) << "The number of tissue types in the mask does not match the number of tissue types in the config file: " << n_tissue << " vs " << param->n_tissue_type;
        return false;
    }
    return true;
}


std::vector<size_t> file_utils::get_size_h5(std::string input_filename, std::string dataset_name)
{
    std::vector<size_t> dims(1,0);
    if(std::filesystem::exists(input_filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "file does not exist: " << input_filename;
        return dims;
    }
    HighFive::File file(input_filename, HighFive::File::ReadOnly);
    if (file.exist(dataset_name) == false)
    {
        BOOST_LOG_TRIVIAL(warning) << "dataset \"" << dataset_name << "\" does not exist in " << input_filename;
        return dims;
    }
    
    HighFive::DataSet dataset = file.getDataSet(dataset_name);
    dims = dataset.getDimensions();
    return dims;    
}


bool file_utils::read_h5(std::string input_filename, void *data, std::string dataset_name, std::string data_type)
{
    if(std::filesystem::exists(input_filename) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "file does not exist: " << input_filename;
        return false;
    }

    HighFive::File file(input_filename, HighFive::File::ReadOnly);
    if (file.exist(dataset_name) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "dataset \"" << dataset_name << "\" does not exist in " << input_filename;
        return false;
    }

    HighFive::DataSet dataset = file.getDataSet(dataset_name);
    if (data_type == "float")
        dataset.read_raw<float>((float*)data);
    else if (data_type == "double")
        dataset.read_raw<double>((double*)data);
    else if (data_type == "uint8_t")
        dataset.read_raw<uint8_t>((uint8_t*)data);
    else
    {
        BOOST_LOG_TRIVIAL(error) << "Data type " << data_type << " is not supported!";
        return false;
    }
    return true;
}


bool file_utils::save_h5(std::string output_filename, void *data, std::vector<size_t> dims, std::string dataset_name, std::string data_type)
{
    std::ostringstream oss;
    for (const auto& elem : dims)
        oss << elem << " ";

    BOOST_LOG_TRIVIAL(info) << "Saving " << dataset_name << " with size = [" << oss.str() << "] to: " << std::filesystem::absolute(output_filename);
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
    // if file exists, open it in read-write mode, otherwise (re-)create it
    auto file_mode = std::filesystem::exists(output_filename) ? HighFive::File::ReadWrite : HighFive::File::Truncate;
    HighFive::File file(output_filename, file_mode);
    // if dataset exists, delete it
    if (file.exist(dataset_name))
        file.unlink(dataset_name);

    if (data_type == "float")
    {
        HighFive::DataSet dataset = file.createDataSet<float>(dataset_name, HighFive::DataSpace(dims));
        dataset.write_raw((float*)data);
    }
    else if (data_type == "double")
    {
        HighFive::DataSet dataset = file.createDataSet<double>(dataset_name, HighFive::DataSpace(dims));
        dataset.write_raw((double*)data);
    }
    else if (data_type == "uint8_t")
    {
        HighFive::DataSet dataset = file.createDataSet<uint8_t>(dataset_name, HighFive::DataSpace(dims));
        dataset.write_raw((uint8_t*)data);
    }
    else
    {
        BOOST_LOG_TRIVIAL(error) << "Data type " << data_type << " is not supported!";
        return false;
    }

    return true;
}


