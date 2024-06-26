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
#define ERR_MSG  "\033[1;31mError:\033[0m "
#define ROUND(x) ((long)((x)+0.5))
#define MAX_RF 256          // maximum number of RF
#define MAX_TE 256          // maximum number of echo times
#define MAX_DEPHASE 256     // maximum number of dephasing
#define MAX_GRADIENT 2048   // maximum number of gradient
#define MAX_TISSUE_TYPE 8   // maximum number of tissue types

typedef struct simulation_parameters
{
    double sample_length[3], scale2grid[3];
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
    uint32_t n_spins, n_timepoints, n_fieldmaps, n_TE, n_RF, n_dephasing, n_gradient, n_sample_length_scales, n_tissue_type;
    uint32_t fieldmap_size[3], seed, max_iterations;
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
        no_gpu(false)
    {
        memset(fieldmap_size,   0, 3*sizeof(fieldmap_size[0])); 
        memset(scale2grid,      0, 3*sizeof(scale2grid[0])); 
        memset(sample_length,   0, 3*sizeof(sample_length[0]));
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
        ss<<"Diffusivity = "; for(int i=0; i<n_tissue_type; i++) ss<<diffusivity[i]<<' '; ss<<"\n";
        ss<<"Cross Tissue Probability =\n"; for(int i=0; i<n_tissue_type; i++) {for(int j=0; j<n_tissue_type; j++) ss<<pXY[j+i*n_tissue_type]<<' '; ss<<'\n';};
        ss<<"RF flip-angle   = "; for(int i=0; i<n_RF; i++) ss<<RF_FA_deg[i]<<' '; ss<<"deg.\n";
        ss<<"RF phase        = "; for(int i=0; i<n_RF; i++) ss<<RF_PH_deg[i]<<' '; ss<<"deg.\n";
        ss<<"RF time         = "; for(int i=0; i<n_RF; i++) ss<<RF_us[i]*timestep_us<<' '; ss<<"us.\n";
        ss<<"dephasing       = "; for(int i=0; i<n_dephasing; i++) ss<<dephasing_deg[i]<<' '; ss<<"deg.\n";
        ss<<"dephasing time  = "; for(int i=0; i<n_dephasing; i++) ss<<dephasing_us[i]*timestep_us<<' '; ss<<"us.\n";
        ss<<"sample length   = "<< sample_length[0] << " x " << sample_length[1] << " x " << sample_length[2] << " m" << '\n';
        ss<<"scale2grid      = "<< scale2grid[0] << " x " << scale2grid[1] << " x " << scale2grid[2] << '\n';
        ss<<"fieldmap size   = "<< fieldmap_size[0] << " x " << fieldmap_size[1] << " x " << fieldmap_size[2] << '\n';
        ss<<"matrix length   = "<< matrix_length << '\n';
        ss<<"dummy scans     = "<< n_dummy_scan<<'\n';
        ss<<"spins           = "<< n_spins<<'\n';
        ss<<"samples scales  = "<< n_sample_length_scales<<'\n';
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
        variables_size_B += sizeof(float) * 3 * spin * n_TE * ((type == "cpu") ? n_sample_length_scales:1);    
        // XYZ0  
        variables_size_B += sizeof(float) * 3 * spin ;     
        // XYZ0_scaled    
        variables_size_B += sizeof(float) * ((type == "gpu") ? 3 * spin : 0);        
        // XYZ1 
        variables_size_B += sizeof(float) * 3 * spin * (enRecordTrajectory ?  n_timepoints * (n_dummy_scan+1) : 1) * ((type == "cpu") ? n_sample_length_scales:1); 
        // T 
        variables_size_B += sizeof(uint8_t) * spin * n_TE * ((type == "cpu") ? n_sample_length_scales:1);
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

        if (mask_exist == false)
        {
            BOOST_LOG_TRIVIAL(error) << "Mask with labeled tissues is not provided";
            return false;
        }
        if (fieldmap_exist == false)
        {
            BOOST_LOG_TRIVIAL(warning) << "Off-resonance map is not provided. Simulation will be performed with a perfect homogeneous field.";
        }
        return true;
    }
} simulation_parameters;



#endif // __SIMULATION_PARAMETERS_H__