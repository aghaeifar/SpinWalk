#include <sstream>
#include <fstream>
#include <filesystem>
#include <iostream>

// boost includes
#include <boost/log/trivial.hpp> 

#include "config_generator.h"
#include "ini.h"

namespace config
{

bool generate_default_config(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output)
{
    BOOST_LOG_TRIVIAL(info) << "Generating default config: " << output << " failed";

    mINI::INIFile file(output);
    mINI::INIStructure ini;

    ini["GENERAL"]["PARENT_CONFIG"] = "";
    ini["GENERAL"]["SEQ_NAME"] = "noname";

    ini["FILES"]["OUTPUT_DIR"] = "./outputs";
    // Off-resonance mapping and masking. The off-resonance map is in Tesla and computed for B0=1T. It will be internally adjusted based on the B0 parameters specified in the \"SIMULATION_PARAMETERS\" section."
    for (size_t i = 0; i < phantoms.size(); i++)
        ini["FILES"]["PHANTOM[" + std::to_string(i) + "]"] = phantoms[i];
    // optional
    ini["FILES"]["XYZ0[0]"] = "";
    ini["FILES"]["XYZ0[1]"] = "";
    ini["FILES"]["M0[0]"] = "";
    ini["FILES"]["M0[1]"] = "";

    // m^2/s
    ini["TISSUE_PARAMETERS"]["DIFFUSIVITY[0]"] = "1.0e-9";
    ini["TISSUE_PARAMETERS"]["DIFFUSIVITY[1]"] = "1.0e-9";
    // Probability to diffuse from tissue X to tissue Y (float). X and Y are taken from the values in the mask
    ini["TISSUE_PARAMETERS"]["P_XY[0]"] = "1.0 0.0";
    ini["TISSUE_PARAMETERS"]["P_XY[1]"] = "0.0 1.0";
    // T1 and T2 in millisecond (float). Negative value to exclude it from the simulation
    ini["TISSUE_PARAMETERS"]["T1[0]"] = "2200";
    ini["TISSUE_PARAMETERS"]["T1[1]"] = "2200";
    ini["TISSUE_PARAMETERS"]["T2[0]"] = "41";
    ini["TISSUE_PARAMETERS"]["T2[1]"] = "41";

    // repetition time in microsecond (integer)
    ini["SCAN_PARAMETERS"]["TR"] = std::to_string(TE_us + timestep_us);
    // echo time in microsecond (integer)
    ini["SCAN_PARAMETERS"]["TE[0]"] = std::to_string(TE_us);
    // RF Flip angle in degree (float)
    ini["SCAN_PARAMETERS"]["RF_FA[0]"] = "90.0";
    // RF Phase in degree (float). Note PHASE_CYCLING will be added to the phase of the first RF
    ini["SCAN_PARAMETERS"]["RF_PH[0]"] = "0.0";
    // Time to apply RF in microsecond (integer). The first RF start time is always 0.0
    ini["SCAN_PARAMETERS"]["RF_T[0]"] = "0";

    // Dephasing in degree (float). The initial spin in the population will experience a dephasing of 0.0 degrees. Dephasing will then progressively increase in a linear manner up to the final spin, which will undergo dephasing as specified by the given parameter
    ini["SCAN_PARAMETERS"]["DEPHASING[0]"] = "";
    // Time to apply dephasing in microsecond (integer).
    ini["SCAN_PARAMETERS"]["DEPHASING_T[0]"] = "";
    // Gradient in mT/m for each axis (float). Each sample is active for one TIME_STEP
    ini["SCAN_PARAMETERS"]["GRADIENT_XYZ[0]"] = "";
    // Time to apply gradient in micro-second (integer).
    ini["SCAN_PARAMETERS"]["GRADIENT_T[0]"] = "";
    // time intervals per random-walk in micro-second (integer)
    ini["SCAN_PARAMETERS"]["TIME_STEP"] = std::to_string(timestep_us);
    // number of dummy scans to reach steady state. The first RF pulse (RF_FA[0]) is used for excitation in dummy scans. If negative, it will be set to 5T1/TR.
    ini["SCAN_PARAMETERS"]["DUMMY_SCAN"] = "0";
    // Phase cycling in degrees
    ini["SCAN_PARAMETERS"]["PHASE_CYCLING"] = "0";

    // static magnetic field in Tesla, set to 0 for no field.
    ini["SIMULATION_PARAMETERS"]["B0"] = "9.4";
    // use 0 for random seed generation, otherwise use a positive integer to make a reproducible simulation
    ini["SIMULATION_PARAMETERS"]["SEED"] = "0";
    ini["SIMULATION_PARAMETERS"]["NUMBER_OF_SPINS"] = "1e5";
    // if 0, spins will not cross volume FoV. If 1, spins which cross FoV, will enter from the other side of the FoV
    ini["SIMULATION_PARAMETERS"]["CROSS_FOV"] = "0";
    // if 1, spins random-walk will be stored in XYZ1 file
    ini["SIMULATION_PARAMETERS"]["RECORD_TRAJECTORY"] = "0";
    // maximum number of iterations that is allowed to generate random-walk. If spin can not move yet (e.g. because of restricted boundaries), it is considered lost and magnetization is set to zero
    ini["SIMULATION_PARAMETERS"]["MAX_ITERATIONS"] = "1e4";
    // SCALE WHAT? 0: FOV, 1: GRADIENT
    ini["SIMULATION_PARAMETERS"]["WHAT_TO_SCALE"] = "0";
    // scale PHANTOM length to simulate different sample sizes
    ini["SIMULATION_PARAMETERS"]["SCALE[0]"] = "1.0";

    
    auto output_file = std::filesystem::absolute(std::filesystem::path(output));
    try {
        std::filesystem::create_directories(output_file.parent_path());
    } catch (const std::exception& e) {
        BOOST_LOG_TRIVIAL(error) << "Creating directory " << output_file.parent_path().string() << " failed. " << e.what();
        return false;
    }

    if(file.generate(ini, true) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Failed to write config file: " << output;
        return false;
    }
    return true;
}


bool generate_gre(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output)
{
    auto output_file = std::filesystem::path(output);
    auto default_config_output = output_file.parent_path() / "default_config.ini";
    if(generate_default_config(TE_us, timestep_us, phantoms, default_config_output.string()) == false)
        return false;

    mINI::INIFile file(output);
    mINI::INIStructure ini;
    
    ini["GENERAL"]["PARENT_CONFIG"] = default_config_output.string();
    ini["GENERAL"]["SEQ_NAME"] = "gre";

    for (size_t i = 0; i < phantoms.size(); i++)
        ini["FILES"]["PHANTOM[" + std::to_string(i) + "]"] = phantoms[i];

    ini["SCAN_PARAMETERS"]["TR"]        = std::to_string(TE_us + timestep_us);
    ini["SCAN_PARAMETERS"]["TE[0]"]     = std::to_string(TE_us);
    ini["SCAN_PARAMETERS"]["RF_FA[0]"]  = "90.0";
    ini["SCAN_PARAMETERS"]["RF_PH[0]"]  = "0";
    ini["SCAN_PARAMETERS"]["RF_T[0]"]   = "0";
    ini["SCAN_PARAMETERS"]["TIME_STEP"] = std::to_string(timestep_us);

    if(file.generate(ini, true) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Failed to write config file: " << output;
        return false;
    }
    return true;
}

bool generate_se(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output)
{
    auto output_file = std::filesystem::path(output);
    auto default_config_output = output_file.parent_path() / "default_config.ini";
    if(generate_default_config(TE_us, timestep_us, phantoms, default_config_output.string()) == false)
        return false;

    mINI::INIFile file(output);
    mINI::INIStructure ini;
    
    ini["GENERAL"]["PARENT_CONFIG"] = default_config_output.string();
    ini["GENERAL"]["SEQ_NAME"] = "se";

    for (size_t i = 0; i < phantoms.size(); i++)
        ini["FILES"]["PHANTOM[" + std::to_string(i) + "]"] = phantoms[i];

    ini["SCAN_PARAMETERS"]["TR"]        = std::to_string(TE_us + timestep_us);
    ini["SCAN_PARAMETERS"]["TE[0]"]     = std::to_string(TE_us);
    ini["SCAN_PARAMETERS"]["RF_FA[0]"]  = "90.0";
    ini["SCAN_PARAMETERS"]["RF_FA[1]"]  = "180.0";
    ini["SCAN_PARAMETERS"]["RF_PH[0]"]  = "0";
    ini["SCAN_PARAMETERS"]["RF_PH[1]"]  = "90";
    ini["SCAN_PARAMETERS"]["RF_T[0]"]   = "0";
    ini["SCAN_PARAMETERS"]["RF_T[1]"]   = std::to_string(TE_us/2);
    ini["SCAN_PARAMETERS"]["TIME_STEP"] = std::to_string(timestep_us);


    if(file.generate(ini, true) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Failed to write config file: " << output;
        return false;
    }
    return true;
}

bool generate_bssfp(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output)
{
    auto output_file = std::filesystem::path(output);
    auto default_config_output = output_file.parent_path() / "default_config.ini";
    if(generate_default_config(TE_us, timestep_us, phantoms, default_config_output.string()) == false)
        return false;

    mINI::INIFile file(output);
    mINI::INIStructure ini;
    
    ini["GENERAL"]["PARENT_CONFIG"] = default_config_output.string();
    ini["GENERAL"]["SEQ_NAME"] = "bssfp";

    for (size_t i = 0; i < phantoms.size(); i++)
        ini["FILES"]["PHANTOM[" + std::to_string(i) + "]"] = phantoms[i];

    ini["SCAN_PARAMETERS"]["TR"]            = std::to_string(TE_us*2);
    ini["SCAN_PARAMETERS"]["TE[0]"]         = std::to_string(TE_us);
    ini["SCAN_PARAMETERS"]["RF_FA[0]"]      = "16.0";
    ini["SCAN_PARAMETERS"]["RF_PH[0]"]      = "0";
    ini["SCAN_PARAMETERS"]["RF_T[0]"]       = "0";
    ini["SCAN_PARAMETERS"]["TIME_STEP"]     = std::to_string(timestep_us);
    ini["SCAN_PARAMETERS"]["DUMMY_SCAN"]    = "-1";
    ini["SCAN_PARAMETERS"]["PHASE_CYCLING"] = "180";


    if(file.generate(ini, true) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "Failed to write config file: " << output;
        return false;
    }
    return true;
}

} // namespace config