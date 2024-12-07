#include <sstream>
#include <fstream>

// boost includes
#include <boost/log/trivial.hpp> 

#include "config_generator.h"
#include "ini.h"

namespace config
{

bool generate_default_config(std::string file_name)
{
    std::stringstream ss;
    ss << "[GENERAL]" << "\n";
    ss << "PARENT_CONFIG = " << "\n";
    ss << "SEQ_NAME = noname" << "\n\n";

    ss << "[FILES]" << "\n";
    ss << "; mandatory" << "\n";
    ss << "OUTPUT_DIR   = ./outputs" << "\n";
    ss << "; Off-resonance mapping and masking. The off-resonance map is in Tesla and computed for B0=1T. It will be internally adjusted based on the B0 parameters specified in the \"SIMULATION_PARAMETERS\" section." << "\n";
    ss << "; The mask defines tissue types. Tissues are labeled from 0 to N-1 for N different tissues. The mask is represented as a 3D matrix with the same dimensions as the Off-resonance map" << "\n";
    ss << "PHANTOM[0]  = ./phantoms/phantom_0.h5" << "\n";
    ss << "; optional" << "\n";
    ss << "XYZ0[0]   = " << "\n";
    ss << "XYZ0[1]   = " << "\n";
    ss << "M0[0]   = " << "\n";
    ss << "M0[1]   = " << "\n\n";

    ss << "[TISSUE_PARAMETERS]" << "\n";
    ss << "; m^2/s (float)" << "\n";
    ss << "DIFFUSIVITY[0] = 1.0e-9" << "\n";
    ss << "DIFFUSIVITY[1] = 1.0e-9" << "\n";
    ss << "; Probability to diffuse from tissue X to tissue Y (float). X and Y are taken from the values in the mask" << "\n";
    ss << "P_XY[0] = 1.0 0.0" << "\n";
    ss << "P_XY[1] = 0.0 1.0" << "\n";
    ss << "; T1 and T2 in millisecond (float). Negative value to exclude it from the simulation" << "\n";
    ss << "T1[0] = 2200" << "\n";
    ss << "T1[1] = 2200" << "\n";
    ss << "T2[0] = 41" << "\n";
    ss << "T2[1] = 41" << "\n\n";

    ss << "[SCAN_PARAMETERS]" << "\n";
    ss << "; repetition time in microsecond (integer)" << "\n";
    ss << "TR = 10e3" << "\n";
    ss << "; echo time in microsecond (integer)" << "\n";
    ss << "TE[0] = 5e3" << "\n";
    ss << "TE[1] = 6e3" << "\n";
    ss << "TE[2] = 7e3" << "\n";
    ss << "; RF Flip angle in degree (float)" << "\n";
    ss << "RF_FA[0] = 15.0" << "\n";
    ss << "RF_FA[1] = 0.0" << "\n";
    ss << "RF_FA[2] = 0.0" << "\n";
    ss << "; RF Phase in degree (float). Note PHASE_CYCLING will be added to the phase of the first RF" << "\n";
    ss << "RF_PH[0] = 0.0" << "\n";
    ss << "RF_PH[1] = 0.0" << "\n";
    ss << "RF_PH[2] = 0.0" << "\n";
    ss << "; Time to apply RF in microsecond (integer). The first RF start time is always 0.0" << "\n";
    ss << "RF_T[0] = 0" << "\n";
    ss << "RF_T[1] = 100e3" << "\n";
    ss << "RF_T[2] = 200e3" << "\n";
    ss << "; Dephasing in degree (float). The initial spin in the population will experience a dephasing of 0.0 degrees. Dephasing will then progressively increase in a linear manner up to the final spin, which will undergo dephasing as specified by the given parameter" << "\n";
    ss << "DEPHASING[0] = " << "\n";
    ss << "DEPHASING[1] = " << "\n";
    ss << "DEPHASING[2] = " << "\n";
    ss << "; Time to apply dephasing in microsecond (integer)." << "\n";
    ss << "DEPHASING_T[0] = " << "\n";
    ss << "DEPHASING_T[1] = " << "\n";
    ss << "DEPHASING_T[2] = " << "\n";
    ss << "; Gradient in mT/m for each axis (float). Each sample is active for one TIME_STEP" << "\n";
    ss << "; GRADIENT_XYZ[x] = 1.0 2.2 1.5" << "\n";
    ss << "GRADIENT_XYZ[0] = " << "\n";
    ss << "GRADIENT_XYZ[1] = " << "\n";
    ss << "GRADIENT_XYZ[2] = " << "\n";
    ss << "; Time to apply gradient in micro-second (integer)." << "\n";
    ss << "GRADIENT_T[0] = " << "\n";
    ss << "GRADIENT_T[1] = " << "\n";
    ss << "GRADIENT_T[2] = " << "\n";
    ss << "; time intervals per random-walk in micro-second (integer)" << "\n";
    ss << "TIME_STEP  = 50" << "\n";
    ss << "; number of dummy scans to reach steady state. The first RF pulse (RF_FA[0]) is used for excitation in dummy scans. If negative, it will be set to 5T1/TR." << "\n";
    ss << "DUMMY_SCAN  = 0" << "\n";
    ss << "; Phase cycling in degrees" << "\n";
    ss << "PHASE_CYCLING = 0" << "\n\n";

    ss << "[SIMULATION_PARAMETERS]" << "\n";
    ss << "; static magnetic field in Tesla, set to 0 for no field. " << "\n";
    ss << "B0   = 9.4" << "\n";
    ss << "; use 0 for random seed generation, otherwise use a positive integer to make a reproducible simulation" << "\n";
    ss << "SEED = 0" << "\n";
    ss << "NUMBER_OF_SPINS = 1e5" << "\n";
    ss << "; if 0, spins will not cross volume FoV. If 1, spins which cross FoV, will enter from the other side of the FoV" << "\n";
    ss << "CROSS_FOV = 0" << "\n";
    ss << "; if 1, spins random-walk will be stored in XYZ1 file" << "\n";
    ss << "RECORD_TRAJECTORY = 0" << "\n";
    ss << "; maximum number of iterations that is allowed to generate random-walk. If spin can not move yet (e.g. because of restricted boundaries), it is considered lost and magnetization is set to zero" << "\n";
    ss << "MAX_ITERATIONS = 1e4" << "\n";
    ss << "; scale PHANTOM length to simulate different sample sizes" << "\n";
    ss << "FOV_SCALE[0] = 1.0" << "\n";
    
    std::ofstream myFile(file_name);
    if (myFile.is_open()) {
        myFile << ss.str();
        myFile.close();
    } else {
       BOOST_LOG_TRIVIAL(error) << "Wriing default config to file " << file_name << " failed";
       return false;
    }

    return true;
}


bool generate_gre(std::string file_name)
{
     return true;
}

bool generate_se(std::string file_name)
{
     return true;
}

bool generate_bssfp(std::string file_name)
{
     return true;
}

bool generate_dwi(std::string file_name)
{
    return true;
}

} // namespace config