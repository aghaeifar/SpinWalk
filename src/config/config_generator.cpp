#include <sstream>
#include <fstream>
#include <filesystem>
#include <iostream>

// boost includes
#include <boost/log/trivial.hpp> 

#include "config_generator.h"

namespace config
{

bool config_generator::generate_default_config(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms)
{
    ini_parent.clear();
    ini_parent["GENERAL"]["PARENT_CONFIG"] = "";
    ini_parent["GENERAL"]["SEQ_NAME"] = "noname";

    ini_parent["FILES"]["OUTPUT_DIR"] = "./outputs";
    // Off-resonance mapping and masking. The off-resonance map is in Tesla and computed for B0=1T. It will be internally adjusted based on the B0 parameters specified in the \"SIMULATION_PARAMETERS\" section."
    for (size_t i = 0; i < phantoms.size(); i++)
        ini_parent["FILES"]["PHANTOM[" + std::to_string(i) + "]"] = phantoms[i];
    // optional
    ini_parent["FILES"]["XYZ0[0]"] = "";
    ini_parent["FILES"]["XYZ0[1]"] = "";
    ini_parent["FILES"]["M0[0]"] = "";
    ini_parent["FILES"]["M0[1]"] = "";

    // m^2/s
    ini_parent["TISSUE_PARAMETERS"]["DIFFUSIVITY[0]"] = "1.0e-9";
    ini_parent["TISSUE_PARAMETERS"]["DIFFUSIVITY[1]"] = "1.0e-9";
    // Probability to diffuse from tissue X to tissue Y (float). X and Y are taken from the values in the mask
    ini_parent["TISSUE_PARAMETERS"]["P_XY[0]"] = "1.0 0.0";
    ini_parent["TISSUE_PARAMETERS"]["P_XY[1]"] = "0.0 1.0";
    // T1 and T2 in millisecond (float). Negative value to exclude it from the simulation
    ini_parent["TISSUE_PARAMETERS"]["T1[0]"] = "2200";
    ini_parent["TISSUE_PARAMETERS"]["T1[1]"] = "2200";
    ini_parent["TISSUE_PARAMETERS"]["T2[0]"] = "41";
    ini_parent["TISSUE_PARAMETERS"]["T2[1]"] = "41";

    // repetition time in microsecond (integer)
    ini_parent["SCAN_PARAMETERS"]["TR"] = std::to_string(TE_us + timestep_us);
    // echo time in microsecond (integer)
    ini_parent["SCAN_PARAMETERS"]["TE[0]"] = std::to_string(TE_us);
    // RF Flip angle in degree (float)
    ini_parent["SCAN_PARAMETERS"]["RF_FA[0]"] = "90.0";
    // RF Phase in degree (float). Note PHASE_CYCLING will be added to the phase of the first RF
    ini_parent["SCAN_PARAMETERS"]["RF_PH[0]"] = "0.0";
    // Time to apply RF in microsecond (integer). The first RF start time is always 0.0
    ini_parent["SCAN_PARAMETERS"]["RF_T[0]"] = "0";

    // Dephasing in degree (float). The initial spin in the population will experience a dephasing of 0.0 degrees. Dephasing will then progressively increase in a linear manner up to the final spin, which will undergo dephasing as specified by the given parameter
    ini_parent["SCAN_PARAMETERS"]["DEPHASING[0]"] = "";
    // Time to apply dephasing in microsecond (integer).
    ini_parent["SCAN_PARAMETERS"]["DEPHASING_T[0]"] = "";
    // Gradient in mT/m for each axis (float). Each sample is active for one TIME_STEP
    ini_parent["SCAN_PARAMETERS"]["GRADIENT_XYZ[0]"] = "";
    // Time to apply gradient in micro-second (integer).
    ini_parent["SCAN_PARAMETERS"]["GRADIENT_T[0]"] = "";
    // time intervals per random-walk in micro-second (integer)
    ini_parent["SCAN_PARAMETERS"]["TIME_STEP"] = std::to_string(timestep_us);
    // number of dummy scans to reach steady state. The first RF pulse (RF_FA[0]) is used for excitation in dummy scans. If negative, it will be set to 5T1/TR.
    ini_parent["SCAN_PARAMETERS"]["DUMMY_SCAN"] = "0";
    // Phase cycling in degrees
    ini_parent["SCAN_PARAMETERS"]["PHASE_CYCLING"] = "0";

    // static magnetic field in Tesla, set to 0 for no field.
    ini_parent["SIMULATION_PARAMETERS"]["B0"] = "9.4";
    // use 0 for random seed generation, otherwise use a positive integer to make a reproducible simulation
    ini_parent["SIMULATION_PARAMETERS"]["SEED"] = "0";
    ini_parent["SIMULATION_PARAMETERS"]["NUMBER_OF_SPINS"] = "1e5";
    // if 0, spins will not cross volume FoV. If 1, spins which cross FoV, will enter from the other side of the FoV
    ini_parent["SIMULATION_PARAMETERS"]["CROSS_FOV"] = "0";
    // if 1, spins random-walk will be stored in XYZ1 file
    ini_parent["SIMULATION_PARAMETERS"]["RECORD_TRAJECTORY"] = "0";
    // maximum number of iterations that is allowed to generate random-walk. If spin can not move yet (e.g. because of restricted boundaries), it is considered lost and magnetization is set to zero
    ini_parent["SIMULATION_PARAMETERS"]["MAX_ITERATIONS"] = "1e4";
    // SCALE WHAT? 0: FOV, 1: GRADIENT
    ini_parent["SIMULATION_PARAMETERS"]["WHAT_TO_SCALE"] = "0";
    // scale PHANTOM length to simulate different sample sizes
    ini_parent["SIMULATION_PARAMETERS"]["SCALE[0]"] = "1.0";

    return true;
}


bool config_generator::generate_gre(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output)
{

    if(generate_default_config(TE_us, timestep_us, phantoms) == false)
        return false;

    ini.clear();
    add_param("GENERAL", "SEQ_NAME", "gre");

    for (size_t i = 0; i < phantoms.size(); i++)
        add_param("FILES", "PHANTOM[" + std::to_string(i) + "]", phantoms[i]);

    add_param("SCAN_PARAMETERS", "TR", std::to_string(TE_us + timestep_us));
    add_param("SCAN_PARAMETERS", "TE[0]", std::to_string(TE_us));
    add_param("SCAN_PARAMETERS", "RF_FA[0]", "90.0");
    add_param("SCAN_PARAMETERS", "RF_PH[0]", "0");
    add_param("SCAN_PARAMETERS", "RF_T[0]", "0");
    add_param("SCAN_PARAMETERS", "TIME_STEP", std::to_string(timestep_us));

    if(output.empty() == false)
        if(write_ini(output) == false)
            return false;

    return true;
}


bool config_generator::generate_se(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output)
{
    if(generate_default_config(TE_us, timestep_us, phantoms) == false)
        return false;

    ini.clear();
    add_param("GENERAL", "SEQ_NAME", "se");

    for (size_t i = 0; i < phantoms.size(); i++)
        add_param("FILES", "PHANTOM[" + std::to_string(i) + "]", phantoms[i]);

    add_param("SCAN_PARAMETERS", "TR", std::to_string(TE_us + timestep_us));
    add_param("SCAN_PARAMETERS", "TE[0]", std::to_string(TE_us));
    add_param("SCAN_PARAMETERS", "RF_FA[0]", "90.0");
    add_param("SCAN_PARAMETERS", "RF_FA[1]", "180.0");
    add_param("SCAN_PARAMETERS", "RF_PH[0]", "0");
    add_param("SCAN_PARAMETERS", "RF_PH[1]", "90");
    add_param("SCAN_PARAMETERS", "RF_T[0]", "0");
    add_param("SCAN_PARAMETERS", "RF_T[1]", std::to_string(TE_us/2));
    add_param("SCAN_PARAMETERS", "TIME_STEP", std::to_string(timestep_us));

    if(output.empty() == false)
        if(write_ini(output) == false)
            return false;

    return true;
}


bool config_generator::generate_bssfp(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output)
{
    if(generate_default_config(TE_us, timestep_us, phantoms) == false)
        return false;

    ini.clear();
    add_param("GENERAL", "SEQ_NAME", "bssfp");

    for (size_t i = 0; i < phantoms.size(); i++)
        add_param("FILES", "PHANTOM[" + std::to_string(i) + "]", phantoms[i]);

    add_param("SCAN_PARAMETERS", "TR", std::to_string(TE_us*2));
    add_param("SCAN_PARAMETERS", "TE[0]", std::to_string(TE_us));
    add_param("SCAN_PARAMETERS", "RF_FA[0]", "16.0");
    add_param("SCAN_PARAMETERS", "RF_PH[0]", "0");
    add_param("SCAN_PARAMETERS", "RF_T[0]", "0");
    add_param("SCAN_PARAMETERS", "TIME_STEP", std::to_string(timestep_us));
    add_param("SCAN_PARAMETERS", "DUMMY_SCAN", "-1");
    add_param("SCAN_PARAMETERS", "PHASE_CYCLING", "180");

    if(output.empty() == false)
        if(write_ini(output) == false)
            return false;
    
    return true;
}


bool config_generator::write_ini(std::string output){
    auto output_file = std::filesystem::absolute(std::filesystem::path(output));
    try {
        std::filesystem::create_directories(output_file.parent_path());
    } catch (const std::exception& e) {
        BOOST_LOG_TRIVIAL(error) << "Creating directory " << output_file.parent_path().string() << " failed. " << e.what();
        return false;
    }
    auto default_config_output = output_file.parent_path() / "default_config.ini";
    add_param("GENERAL", "PARENT_CONFIG", default_config_output.string().c_str());

    mINI::INIFile file(output_file.string());
    if(file.generate(ini, true) == false){
        BOOST_LOG_TRIVIAL(error) << "Failed to write config file: " << output;
        return false;
    }

    mINI::INIFile file_parent(default_config_output.string());
    if(file_parent.generate(ini_parent, true) == false){
        BOOST_LOG_TRIVIAL(error) << "Failed to write config file: " << default_config_output;
        return false;
    }

    return true;
}

void config_generator::add_param(std::string section, std::string key, std::string value){
    ini[section][key] = value;
}

} // namespace config