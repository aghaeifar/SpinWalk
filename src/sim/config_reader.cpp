#include <filesystem>
#include <iostream>

#include "../definitions.h"
#include "config_reader.h"
#include "simulation_parameters.cuh"
#include "ini.h"

// boost includes
#include <boost/log/trivial.hpp> 

namespace sim{

bool config_reader::prepare(std::string config_filename, simulation_parameters *param)
{
    cleanup();
    this->param = param;
    this->config_filename = config_filename;
    if (this->read(config_filename) == false)
        return false;
    if (this->check() == false)
        return false;
    return true;
}

void config_reader::cleanup()
{
    seq_name = "";
    output_dir = "";
    config_filename = "";
    fov_scales.clear();
    output_files.clear();
    for (auto &files : files_container )
        files.second.clear();
    param = nullptr;
}

bool config_reader::read(std::string config_filename_path)
{
    std::stringstream ss;
    if (std::filesystem::exists(config_filename_path) == false){
        BOOST_LOG_TRIVIAL(error) << "Config-file does not exist: " << config_filename_path;
        return false;
    }
    
    std::string filename = std::filesystem::canonical(std::filesystem::path(config_filename_path)).filename().string();
    BOOST_LOG_TRIVIAL(info) << "Reading config: " << config_filename_path;
    mINI::INIFile file(config_filename_path);
    mINI::INIStructure ini;
    if(file.read(ini) == false){
        std::cout << ERR_MSG << "Failed to read config file: " << config_filename_path << "\n";
        return false;
    }
    
    // ============== reading section GENERAL ==============
    if(ini["GENERAL"]["PARENT_CONFIG"].empty() == false) {
        std::filesystem::path parent_config(ini["GENERAL"]["PARENT_CONFIG"]);   
        if (parent_config.is_relative()) // if parent_config is relative, make it absolute
            parent_config = std::filesystem::absolute(config_filename_path).parent_path() / parent_config;
        if (this->read(parent_config.string()) == false)
            return false;
        BOOST_LOG_TRIVIAL(info) << "Back to reading config: " << config_filename_path;
    }

    seq_name = ini["GENERAL"]["SEQ_NAME"].empty() ? seq_name : ini["GENERAL"]["SEQ_NAME"];
    
    // ============== reading section FILES ==============
    for (const auto& str : {std::string("PHANTOM"), std::string("XYZ0"), std::string("M0")}) 
    {
        if(ini["FILES"].has(str + "[0]"))
            files_container.at(str).clear();  
        for (uint16_t i = 0; ini["FILES"][str + "[" + std::to_string(i) + "]"].empty() == false ; i++) 
            files_container.at(str).push_back(ini["FILES"][str + "[" + std::to_string(i) + "]"]);  
        // make paths absolute
        for (auto& path : files_container.at(str))
            if (std::filesystem::path(path).is_relative()) 
                path = std::filesystem::canonical(std::filesystem::absolute(config_filename_path).parent_path() / path).string();
    }
    
    output_dir = ini["FILES"]["OUTPUT_DIR"].empty() ? output_dir : ini["FILES"]["OUTPUT_DIR"]; 
    if (std::filesystem::path(output_dir).is_relative())
        output_dir = std::filesystem::canonical(std::filesystem::absolute(config_filename).parent_path() / output_dir).string();
     
    // ============== reading section SCAN_PARAMETERS ==============
    param->TR_us = ini["SCAN_PARAMETERS"]["TR"].empty() ? param->TR_us : std::stoi(ini["SCAN_PARAMETERS"]["TR"]);
    param->timestep_us = ini["SCAN_PARAMETERS"]["TIME_STEP"].empty() ? param->timestep_us : std::stoi(ini["SCAN_PARAMETERS"]["TIME_STEP"]);

    // ---------------- Echo times ----------------      
    std::vector<int32_t> TE_us; 
    for(uint16_t i=0; ini["SCAN_PARAMETERS"]["TE[" + std::to_string(i) + "]"].empty() == false ; i++) 
        TE_us.push_back(std::stoi(ini["SCAN_PARAMETERS"]["TE[" + std::to_string(i) + "]"]) / param->timestep_us );    
    if(TE_us.size() > 0)
        param->TE_us.data = TE_us;
   
    // ---------------- RF pulses (start times, flip angles, phases ) ----------------
    // RF start times
    std::vector<int32_t> RF_us; 
    for(uint16_t i=0; ini["SCAN_PARAMETERS"]["RF_T[" + std::to_string(i) + "]"].empty() == false ; i++) 
        RF_us.push_back(std::stoi(ini["SCAN_PARAMETERS"]["RF_T[" + std::to_string(i) + "]"]) / param->timestep_us );       
    if (RF_us.size() > 0)
        param->RF_us.data = RF_us;  
    
    // RF flip angles
    std::vector<float> RF_FA; 
    for(uint16_t i=0; ini["SCAN_PARAMETERS"]["RF_FA[" + std::to_string(i) + "]"].empty() == false ; i++) 
        RF_FA.push_back(std::stof(ini["SCAN_PARAMETERS"]["RF_FA[" + std::to_string(i) + "]"])); 
    if (RF_FA.size() > 0)
        param->RF_FA_deg.data = RF_FA;

    // RF phases
    std::vector<float> RF_PH; 
    for(uint16_t i=0; ini["SCAN_PARAMETERS"]["RF_PH[" + std::to_string(i) + "]"].empty() == false ; i++) 
        RF_PH.push_back(std::stof(ini["SCAN_PARAMETERS"]["RF_PH[" + std::to_string(i) + "]"])); 
    if (RF_PH.size() > 0)
        param->RF_PH_deg.data = RF_PH;

    // ---------------- dephasing (start times, Flip angles ) ----------------
    // Dephase start times
    std::vector<int32_t> dephasing_us; 
    for(uint16_t i=0; ini["SCAN_PARAMETERS"]["DEPHASING_T[" + std::to_string(i) + "]"].empty() == false ; i++) 
        dephasing_us.push_back(std::stoi(ini["SCAN_PARAMETERS"]["DEPHASING_T[" + std::to_string(i) + "]"]) / param->timestep_us );       
    if (dephasing_us.size() > 0)
        param->dephasing_us.data = dephasing_us; 
     
    // Dephase flip angles
    std::vector<float> dephasing_deg; 
    for(uint16_t i=0; ini["SCAN_PARAMETERS"]["DEPHASING[" + std::to_string(i) + "]"].empty() == false ; i++) 
        dephasing_deg.push_back(std::stof(ini["SCAN_PARAMETERS"]["DEPHASING[" + std::to_string(i) + "]"])); 
    if (dephasing_deg.size() > 0)
        param->dephasing_deg.data = dephasing_deg;
     
    // ---------------- Gradients (start times, strength (T/m) ) ----------------
    // Gradient start times
    std::vector<int32_t> gradient_us; 
    for(uint16_t i=0; ini["SCAN_PARAMETERS"]["GRADIENT_T[" + std::to_string(i) + "]"].empty() == false ; i++) 
        gradient_us.push_back(std::stoi(ini["SCAN_PARAMETERS"]["GRADIENT_T[" + std::to_string(i) + "]"]) / param->timestep_us );       
    if (gradient_us.size() > 0)
        param->gradient_us.data = gradient_us; 
     
    // Gradient strength
    std::vector<float> gradient_mTm; 
    for(uint16_t i=0; ini["SCAN_PARAMETERS"]["GRADIENT_XYZ[" + std::to_string(i) + "]"].empty() == false ; i++)     {
        std::istringstream iss(ini["SCAN_PARAMETERS"]["GRADIENT_XYZ[" + std::to_string(i) + "]"]);
        std::vector<double> values{std::istream_iterator<double>(iss), std::istream_iterator<double>()};
        gradient_mTm.insert(gradient_mTm.end(), values.begin(), values.end());
    }
    if (gradient_mTm.size() > 0)
        param->gradient_mTm.data = gradient_mTm;
     
    // ============== reading section SCAN_PARAMETERS ==============
    param->n_dummy_scan = ini["SCAN_PARAMETERS"]["DUMMY_SCAN"].empty() ? param->n_dummy_scan : std::stoi(ini["SCAN_PARAMETERS"]["DUMMY_SCAN"]);
    param->phase_cycling = ini["SCAN_PARAMETERS"]["PHASE_CYCLING"].empty() ? param->phase_cycling : std::stof(ini["SCAN_PARAMETERS"]["PHASE_CYCLING"]);
     
    // ============== reading section SIMULATION_PARAMETERS ==============
    param->B0                    = ini["SIMULATION_PARAMETERS"]["B0"].empty() ? param->B0 : std::stof(ini["SIMULATION_PARAMETERS"]["B0"]);
    param->seed                  = ini["SIMULATION_PARAMETERS"]["SEED"].empty() ? param->seed : std::stoi(ini["SIMULATION_PARAMETERS"]["SEED"]);
    param->n_spins               = ini["SIMULATION_PARAMETERS"]["NUMBER_OF_SPINS"].empty() ? param->n_spins : std::stod(ini["SIMULATION_PARAMETERS"]["NUMBER_OF_SPINS"]); // template type must be double since input can be of form scientific notation
    param->enCrossFOV            = ini["SIMULATION_PARAMETERS"]["CROSS_FOV"].empty() ? param->enCrossFOV : std::stoi(ini["SIMULATION_PARAMETERS"]["CROSS_FOV"]);
    param->enRecordTrajectory    = ini["SIMULATION_PARAMETERS"]["RECORD_TRAJECTORY"].empty() ? param->enRecordTrajectory : std::stoi(ini["SIMULATION_PARAMETERS"]["RECORD_TRAJECTORY"]);
    param->max_iterations        = ini["SIMULATION_PARAMETERS"]["MAX_ITERATIONS"].empty() ? param->max_iterations : std::stod(ini["SIMULATION_PARAMETERS"]["MAX_ITERATIONS"]);
     
    // Field of view Scales
    if(ini["SIMULATION_PARAMETERS"].has("FOV_SCALE[0]"))
        fov_scales.clear();
    for (uint16_t i = 0; ini["SIMULATION_PARAMETERS"]["FOV_SCALE[" + std::to_string(i) + "]"].empty() == false ; i++) 
        fov_scales.push_back(std::stod(ini["SIMULATION_PARAMETERS"]["FOV_SCALE[" + std::to_string(i) + "]"]));
    param->n_fov_scale = fov_scales.size(); 
      
    // ============== reading section TISSUE_PARAMETERS ==============
    // Diffusivity
    std::vector<double> diffusivity; 
    for(uint16_t i=0; ini["TISSUE_PARAMETERS"]["DIFFUSIVITY[" + std::to_string(i) + "]"].empty() == false ; i++) 
        diffusivity.push_back(std::stof(ini["TISSUE_PARAMETERS"]["DIFFUSIVITY[" + std::to_string(i) + "]"])); 
    if (diffusivity.size() > 0)
        param->diffusivity.data = diffusivity;
    param->n_substrate = param->diffusivity.data.size();
     
    // T1 & T2
    std::vector<float> T1_ms; 
    for(uint16_t i=0; ini["TISSUE_PARAMETERS"]["T1[" + std::to_string(i) + "]"].empty() == false ; i++) 
        T1_ms.push_back(std::stof(ini["TISSUE_PARAMETERS"]["T1[" + std::to_string(i) + "]"])); 
    if (T1_ms.size() > 0)
        param->T1_ms.data = T1_ms;
     
    std::vector<float> T2_ms; 
    for(uint16_t i=0; ini["TISSUE_PARAMETERS"]["T2[" + std::to_string(i) + "]"].empty() == false ; i++) 
        T2_ms.push_back(std::stof(ini["TISSUE_PARAMETERS"]["T2[" + std::to_string(i) + "]"])); 
    if (T2_ms.size() > 0)
        param->T2_ms.data = T2_ms;
     
    // Cross Tissue Probability
    std::vector<float> pXY; 
    for(uint16_t i=0; ini["TISSUE_PARAMETERS"]["P_XY[" + std::to_string(i) + "]"].empty() == false ; i++) {
        std::istringstream iss(ini["TISSUE_PARAMETERS"]["P_XY[" + std::to_string(i) + "]"]);
        std::vector<double> values{std::istream_iterator<double>(iss), std::istream_iterator<double>()};
        pXY.insert(pXY.end(), values.begin(), values.end());
    }
    if (pXY.size() > 0)
        param->pXY.data = pXY;
         
    // ============== End ==============
    return true;
}

bool config_reader::check()
{
    BOOST_LOG_TRIVIAL(info) << "Checking consistentcy of parameters in config file...";
    // ============== check files ==============
    for ( const auto &files : files_container )
        for (const auto& path : files.second)
            if (std::filesystem::exists(path) == false) {
                BOOST_LOG_TRIVIAL(error) << "File does not exist: " << path;
                return false;
            }

    // resize for consistency
    files_container.at("XYZ0").resize(files_container.at("PHANTOM").size(), "");
    files_container.at("M0").resize(files_container.at("PHANTOM").size(), "");
    
    // generate names for output
    output_files.clear();
    for (const auto& path : files_container.at("PHANTOM")){
        auto f = std::filesystem::path(output_dir) / (seq_name + "_" + std::filesystem::path(path).filename().string());
        output_files.push_back(f.replace_extension(".h5").string());
    }

    // ============== check timings ==============
    // TE
    if (std::is_sorted(param->TE_us.data.begin(), param->TE_us.data.end()) == false || 
        std::adjacent_find(param->TE_us.data.begin(), param->TE_us.data.end()) != param->TE_us.data.end() || 
        param->TE_us.data[0] < 0 || param->TE_us.data.size() == 0)
    {
        std::stringstream ss; std::copy(param->TE_us.data.begin(), param->TE_us.data.end(), std::ostream_iterator<int>(ss, " "));
        BOOST_LOG_TRIVIAL(error) << "TE must exists and be in ascending order and must not have duplicates or negative values: " << ss.str();
        return false;
    }     

    // RF
    if (std::is_sorted(param->RF_us.data.begin(), param->RF_us.data.end()) == false || 
        std::adjacent_find(param->RF_us.data.begin(), param->RF_us.data.end()) != param->RF_us.data.end() || 
        param->RF_us.data[0] != 0 || param->RF_us.data.size() == 0)
    {
        std::stringstream ss; std::copy(param->RF_us.data.begin(), param->RF_us.data.end(), std::ostream_iterator<int>(ss, " "));
        BOOST_LOG_TRIVIAL(error) << "RF times must be in ascending order, starts with 0 and must not have duplicates values: " << ss.str();
        return false;
    }  

    // dephasing   
    if (std::is_sorted(param->dephasing_us.data.begin(), param->dephasing_us.data.end()) == false || 
        std::adjacent_find(param->dephasing_us.data.begin(), param->dephasing_us.data.end()) != param->dephasing_us.data.end())
    {
        std::stringstream ss; std::copy(param->dephasing_us.data.begin(), param->dephasing_us.data.end(), std::ostream_iterator<int>(ss, " "));
        BOOST_LOG_TRIVIAL(error) << "Dephasing Times must be in ascending order and must not have duplicates values: " << ss.str();
        return false;
    }   

    // gradients   
    if (std::is_sorted(param->gradient_us.data.begin(), param->gradient_us.data.end()) == false || 
        std::adjacent_find(param->gradient_us.data.begin(), param->gradient_us.data.end()) != param->gradient_us.data.end())
    {
        std::stringstream ss; std::copy(param->gradient_us.data.begin(), param->gradient_us.data.end(), std::ostream_iterator<int>(ss, " "));
        BOOST_LOG_TRIVIAL(error) << "Gradient times must be in a strickly ascending order and must not have duplicates values: " << ss.str();
        return false;
    }   

    // ============== check size ==============
    if( param->RF_FA_deg.data.size() != param->RF_us.data.size() || param->RF_FA_deg.data.size() != param->RF_PH_deg.data.size())
    {
        BOOST_LOG_TRIVIAL(error) << "RF_FA, RF_PH and RF_us must have the same number of elements " << param->RF_FA_deg.data.size() << " vs " << param->RF_PH_deg.data.size() << " vs " << param->RF_us.data.size();
        return false;
    }

    if(param->dephasing_us.data.size() != param->dephasing_deg.data.size())
    {
        BOOST_LOG_TRIVIAL(error) << "DEPHASING and DEPHASING_T must have the same number of elements " << param->dephasing_deg.data.size() << " vs " << param->dephasing_us.data.size();
        return false;
    }

    if(param->gradient_mTm.data.size() / 3 != param->gradient_us.data.size())
    {
        BOOST_LOG_TRIVIAL(error) << "GRADIENT_XYZ and GRADIENT_T must have the same number of elements " << param->gradient_mTm.data.size() << " vs " << param->gradient_us.data.size();
        return false;
    }

    if (param->T1_ms.data.size() != param->T2_ms.data.size())
    {
        BOOST_LOG_TRIVIAL(error) << "T1 and T2 must have the same number of elements " << param->T1_ms.data.size() << " vs " << param->T2_ms.data.size();
        return false;
    }

    if (param->T1_ms.data.size() != param->diffusivity.data.size())
    {
        BOOST_LOG_TRIVIAL(error) << "T1 and diffusivity must have the same number of elements " << param->T1_ms.data.size() << " vs " << param->diffusivity.data.size();
        return false;
    }

    if (param->T1_ms.data.size() * param->T1_ms.data.size() != param->pXY.data.size())
    {
        BOOST_LOG_TRIVIAL(error) << "T1 and P_XY must have the same number of elements " << param->T1_ms.data.size() << " vs " << param->pXY.data.size();
        return false;
    }

    if (fov_scales.size() == 0)
    {
        BOOST_LOG_TRIVIAL(warning) << "FOV_SCALE is not set! Using default value 1.0";
        fov_scales.push_back(1.0);
    }

    return true;
}

void config_reader::dump() const
{

    mINI::INIFile file(config_filename);
    mINI::INIStructure ini;
    if(file.read(ini) == false){
        std::cout << ERR_MSG << "Failed to read config file: " << config_filename << "\n";
        return;
    }

    for (auto const& it : ini){
        auto const& section = it.first;
        auto const& collection = it.second;
        std::cout << "[" << section << "]" << "\n";
        for (auto const& it2 : collection)
        {
            auto const& key = it2.first;
            auto const& value = it2.second;
            std::cout << key << "=" << value << "\n";
        }
        std::cout << "\n";
    }
}

} // namespace sim