/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: miscellaneous.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#ifndef __SIMULATION_PARAMETERS_H__
#define __SIMULATION_PARAMETERS_H__

#include <random>
#include <algorithm>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <cctype>
#include <cmath>
#include <map>
#include <boost/log/trivial.hpp> 

#define DEG2RAD 0.0174532925199433 // = M_PI/180 
#define RAD2DEG 57.2957795130823

#define B2MB 1048576
#define ROUND(x) ((long)((x)+0.5))
#define MAX_RF 256          // maximum number of RF
#define MAX_TE 256          // maximum number of echo times
#define MAX_DEPHASE 256     // maximum number of dephasing
#define MAX_GRADIENT 2048   // maximum number of gradient
#define MAX_TISSUE_TYPE 8   // maximum number of tissue types

typedef struct simulation_parameters
{
    double fov[3], scale2grid[3];
    double diffusivity[MAX_TISSUE_TYPE];
    float B0, c, s;
    float RF_FA_deg[MAX_RF], RF_PH_deg[MAX_RF];
    float dephasing_deg[MAX_DEPHASE];          
    float gradient_mTm[3*MAX_GRADIENT];     
    float pXY[MAX_TISSUE_TYPE*MAX_TISSUE_TYPE];
    float phase_cycling;
    float T1_ms[MAX_TISSUE_TYPE], T2_ms[MAX_TISSUE_TYPE];
    int32_t timestep_us, TR_us, TE_us[MAX_TE], RF_us[MAX_RF], dephasing_us[MAX_DEPHASE], gradient_us[MAX_GRADIENT];
    int32_t n_dummy_scan;
    uint32_t n_spins, n_timepoints, n_fieldmaps, n_TE, n_RF, n_dephasing, n_gradient, n_fov_scale, n_tissue_type;
    size_t fieldmap_size[3], seed, max_iterations;
    int64_t matrix_length;
    bool enDebug, enCrossFOV, enRecordTrajectory, enProfiling;
    bool fieldmap_exist, mask_exist, no_gpu;
    
    simulation_parameters():
        TR_us(40e3),
        timestep_us(20),
        B0(9.4),
        n_TE(0),
        n_RF(0),
        n_dephasing(0),
        n_gradient(0),
        n_dummy_scan(0),
        n_tissue_type(0),
        max_iterations(9999),
        phase_cycling(0.),
        enDebug(false),
        enCrossFOV(true),
        enRecordTrajectory(false),
        fieldmap_exist(true),
        mask_exist(true),
        enProfiling(false),
        matrix_length(0),
        no_gpu(false)
    {
        memset(fieldmap_size,   0, 3*sizeof(fieldmap_size[0])); 
        memset(scale2grid,      0, 3*sizeof(scale2grid[0])); 
        memset(fov,   0, 3*sizeof(fov[0]));
        memset(TE_us,           0, MAX_TE*sizeof(TE_us[0]));        
        memset(RF_FA_deg,       0, MAX_RF*sizeof(RF_FA_deg[0]));
        memset(RF_us,           0, MAX_RF*sizeof(RF_us[0]));
        memset(RF_PH_deg,       0, MAX_RF*sizeof(RF_PH_deg[0]));
        memset(dephasing_deg,   0, MAX_DEPHASE*sizeof(dephasing_deg[0]));
        memset(dephasing_us,    0, MAX_DEPHASE*sizeof(dephasing_us[0]));
        memset(gradient_mTm,    0, 3*MAX_GRADIENT*sizeof(gradient_mTm[0]));
        memset(gradient_us,     0, MAX_GRADIENT*sizeof(gradient_us[0]));

        std::fill(pXY, pXY + MAX_TISSUE_TYPE*MAX_TISSUE_TYPE, 1.f);
        std::fill(diffusivity, diffusivity + MAX_TISSUE_TYPE, 1.0);
        std::fill(T1_ms, T1_ms + MAX_TISSUE_TYPE, 2000);
        std::fill(T2_ms, T2_ms + MAX_TISSUE_TYPE, 45);
    }

    std::string dump()
    {
        std::stringstream ss;
        ss<<"Use GPU = "<<(no_gpu?"No":"Yes")<<'\n';
        ss<<"B0 = "<<B0<<" T\n";
        ss<<"timestep = "<<timestep_us<<" us.\n";
        ss<<"TR = "<<TR_us/1000.<<" ms.\n";
        ss<<"TE = "; for(int i=0; i<n_TE; i++) ss<<TE_us[i]*timestep_us/1000.<<' '; ss<<"ms.\n";
        ss<<"T2 = "; for(int i=0; i<n_tissue_type; i++) ss<<T2_ms[i]<<' '; ss<<"ms.\n";
        ss<<"T1 = "; for(int i=0; i<n_tissue_type; i++) ss<<T1_ms[i]<<' '; ss<<"ms.\n";
        ss<<"sqrt(2*diffusivity*timestep) = "; for(int i=0; i<n_tissue_type; i++) ss<<diffusivity[i]<<' '; ss<<"\n";
        ss<<"Cross Tissue Probability =\n"; for(int i=0; i<n_tissue_type; i++) {for(int j=0; j<n_tissue_type; j++) ss<<pXY[j+i*n_tissue_type]<<' '; ss<<'\n';};
        ss<<"RF flip-angle   = "; for(int i=0; i<n_RF; i++) ss<<RF_FA_deg[i]<<' '; ss<<"deg.\n";
        ss<<"RF phase        = "; for(int i=0; i<n_RF; i++) ss<<RF_PH_deg[i]<<' '; ss<<"deg.\n";
        ss<<"RF time         = "; for(int i=0; i<n_RF; i++) ss<<RF_us[i]*timestep_us<<' '; ss<<"us.\n";
        ss<<"dephasing       = "; for(int i=0; i<n_dephasing; i++) ss<<dephasing_deg[i]<<' '; ss<<"deg.\n";
        ss<<"dephasing time  = "; for(int i=0; i<n_dephasing; i++) ss<<dephasing_us[i]*timestep_us<<' '; ss<<"us.\n";
        ss<<"FoV             = "<< fov[0] << " x " << fov[1] << " x " << fov[2] << " m" << '\n';
        ss<<"scale2grid      = "<< scale2grid[0] << " x " << scale2grid[1] << " x " << scale2grid[2] << '\n';
        ss<<"fieldmap size   = "<< fieldmap_size[0] << " x " << fieldmap_size[1] << " x " << fieldmap_size[2] << '\n';
        ss<<"matrix length   = "<< matrix_length << '\n';
        ss<<"dummy scans     = "<< n_dummy_scan<<'\n';
        ss<<"spins           = "<< n_spins<<'\n';
        ss<<"FoV scales      = "<< n_fov_scale<<'\n';
        ss<<"timepoints      = "<< n_timepoints<<'\n';
        ss<<"fieldmaps       = "<< n_fieldmaps<<'\n';
        ss<<"max iterations  = "<< max_iterations<<'\n';
        ss<<"Tissue Types    = "<< n_tissue_type << '\n';
        ss<<"Pass FoV        = "<< enCrossFOV << '\n';
        ss<<"Phase cycling   = "<< phase_cycling<<'\n';
        ss<<"Seed            = "<< seed<<'\n';
        ss<<"mask exists     = "<< mask_exist<<'\n';
        ss<<"off-resonance exists   = "<< fieldmap_exist<<'\n';
        ss<<"gradient (x,y,z) mT/m  =\n"; for(int i=0; i<n_gradient; i++) ss<<gradient_mTm[3*i+0]<<' '<<gradient_mTm[3*i+1]<<' '<<gradient_mTm[3*i+2]<<'\n';
        ss<<"gradient  time         =\n"; for(int i=0; i<n_gradient; i++) ss<<gradient_us[i]*timestep_us<<' '; ss<<"us.\n";
        ss<<"Record Trajectory      = "<< enRecordTrajectory << '\n';
        ss<<"Required CPU memory    = "<<get_required_memory(1, "cpu")<<" MB\n";
        ss<<"Required GPU memory    = "<<get_required_memory(1, "gpu")<<" MB\n";

        return ss.str();
    }

    size_t get_required_memory(uint8_t n_device=1, std::string type="gpu")
    {
        std::transform(type.begin(), type.end(), type.begin(), [](unsigned char c){ return std::tolower(c); }); 
        size_t spin = (type == "gpu") ? n_spins/n_device : n_spins;
        // fieldmap and mask
        size_t data_size_MB = 0;
        // off-resonance map
        data_size_MB += fieldmap_exist ? (matrix_length * sizeof(float) / B2MB) : 0;
        // mask
        data_size_MB += mask_exist ? (matrix_length * sizeof(uint8_t) / B2MB) : 0;
        // variables (M0, XYZ0, XYZ0_scaled, XYZ1, M1)
        size_t variables_size_B = 0;
        // M0
        variables_size_B += sizeof(float) * 3 * spin ;     
        // M1     
        variables_size_B += sizeof(float) * 3 * spin * n_TE * ((type == "cpu") ? n_fov_scale:1);    
        // XYZ0  
        variables_size_B += sizeof(float) * 3 * spin ;     
        // XYZ0_scaled    
        variables_size_B += sizeof(float) * ((type == "gpu") ? 3 * spin : 0);        
        // XYZ1 
        variables_size_B += sizeof(float) * 3 * spin * (enRecordTrajectory ?  n_timepoints * (n_dummy_scan+1) : 1) * ((type == "cpu") ? n_fov_scale:1); 
        // T 
        variables_size_B += sizeof(uint8_t) * spin * n_TE * ((type == "cpu") ? n_fov_scale:1);
        size_t variables_size_MB = variables_size_B / B2MB;

        return data_size_MB + variables_size_MB;
    }

    bool prepare()
    {
        c = cosf(RF_FA_deg[0] * DEG2RAD); 
        s = sinf(RF_FA_deg[0] * DEG2RAD); 
        matrix_length = fieldmap_size[0] * fieldmap_size[1] * fieldmap_size[2];
        n_timepoints = TR_us / timestep_us;
        TE_us[n_TE++] = n_timepoints - 1; // add the last timepoint

        if (seed == 0)
            seed = std::random_device{}();

        for (int i = 0; i < n_tissue_type; i++)
            diffusivity[i] = 1e-3 * sqrt(2. * diffusivity[i] * timestep_us); // 1e-3 is the conversion factor from us to seconds for timestep_us

        if (n_dummy_scan < 0)
            n_dummy_scan = 5.0 * T1_ms[0] / TR_us * 1e3;

        if (n_TE == 1)
        {
            BOOST_LOG_TRIVIAL(error) << "Echo time is not set";
            return false;
        }        
        if (mask_exist == false)
        {
            BOOST_LOG_TRIVIAL(error) << "Mask with labeled tissues is not provided";
            return false;
        }
        if (fieldmap_exist == false)
        {
            BOOST_LOG_TRIVIAL(warning) << "Off-resonance map is not provided. Simulation will be performed with a perfect homogeneous field.";
        }
        if(n_fov_scale == 0)
        {
            BOOST_LOG_TRIVIAL(error) << "No FoV scale is provided";
            return false;
        }
        return true;
    }
} simulation_parameters;


// --------------------------------------------------------------------------------------------
// generate default configuration file
// --------------------------------------------------------------------------------------------
inline bool generate_default_config(std::string file_name)
{
    std::stringstream ss;
    ss << "[PARENT]" << "\n";
    ss << "PARENT_CONFIG = " << "\n\n";

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

#endif // __SIMULATION_PARAMETERS_H__