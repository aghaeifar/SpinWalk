/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: miscellaneous.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#ifndef __MISCELLANEOUS_H__
#define __MISCELLANEOUS_H__


#include <algorithm>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <cctype>
#include <cmath>

#define SPINWALK_VERSION_MAJOR 1
#define SPINWALK_VERSION_MINOR 7
#define SPINWALK_VERSION_PATCH 1

#define DEG2RAD 0.0174532925199433 // = M_PI/180 
#define RAD2DEG 57.2957795130823

#define B2MB 1048576
#define ERR_MSG  "\033[1;31mError:\033[0m "
#define ROUND(x) ((long)((x)+0.5))
#define MAX_RF 256          // maximum number of RF
#define MAX_TE 256          // maximum number of echo times
#define MAX_T12 256         // maximum number of relaxation times
#define MAX_DEPHASE 256     // maximum number of dephasing
#define MAX_GRADIENT 256    // maximum number of gradient
#define MAX_TISSUE_TYPE 8   // maximum number of tissue types

typedef struct simulation_parameters
{
    float B0, c, s, c2, s2;
    float T1[MAX_T12], T2[MAX_T12];
    float RF_FA[MAX_RF], RF_PH[MAX_RF]; // refocusing FA
    float dephasing[MAX_DEPHASE]; // dephasing in degree
    float gradient_xyz[3*MAX_GRADIENT]; // gradient in T/m
    float pXY[MAX_TISSUE_TYPE*MAX_TISSUE_TYPE];
    uint32_t RF_ST[MAX_RF], TE[MAX_TE], dephasing_T[MAX_DEPHASE], gradient_T[MAX_GRADIENT]; // refocusing time in dt, echo times in dt, dephasing time in dt
    double sample_length[3], scale2grid[3], diffusion_const, TR, dt;
    float phase_cycling;
    uint32_t n_timepoints, n_sample_length_scales, n_fieldmaps, n_TE, n_RF, n_dephasing, n_gradient, n_T12;
    int32_t n_dummy_scan ;
    uint32_t n_tissue_type;
    uint32_t n_spins, fieldmap_size[3], seed, max_iterations;
    int64_t matrix_length;
    bool enDebug, enCrossBoundry, enRecordTrajectory;
    simulation_parameters():
        TR(0.04),
        dt(5e-5),
        B0(9.4),
        n_TE(0),
        n_RF(0),
        n_T12(0),
        n_dephasing(0),
        n_gradient(0),
        n_dummy_scan(0),
        n_tissue_type(0),
        max_iterations(9999),
        diffusion_const(1.2e-9),
        phase_cycling(0.),
        enDebug(false),
        enCrossBoundry(true),
        enRecordTrajectory(false)
    {
        memset(fieldmap_size, 0, 3*sizeof(fieldmap_size[0])); 
        memset(sample_length, 0, 3*sizeof(sample_length[0]));
        memset(TE, 0, MAX_TE*sizeof(TE[0]));        
        memset(RF_FA, 0, MAX_RF*sizeof(RF_FA[0]));
        memset(RF_ST, 0, MAX_RF*sizeof(RF_ST[0]));
        memset(RF_PH, 0, MAX_RF*sizeof(RF_PH[0]));
        memset(dephasing, 0, MAX_DEPHASE*sizeof(dephasing[0]));
        memset(dephasing_T, 0, MAX_DEPHASE*sizeof(dephasing_T[0]));
        memset(gradient_xyz, 0, 3*MAX_GRADIENT*sizeof(gradient_xyz[0]));
        memset(gradient_T, 0, MAX_GRADIENT*sizeof(gradient_T[0]));
        memset(pXY, 0, MAX_TISSUE_TYPE*MAX_TISSUE_TYPE*sizeof(pXY[0]));

        std::fill(T1, T1 + MAX_T12, 2.2);
        std::fill(T2, T2 + MAX_T12, 0.04);
    }

    std::string dump()
    {
        std::stringstream ss;
        ss<<"B0 = "<<B0<<'\n'<<"dt = "<<dt<<'\n';
        ss<<"TR = "<<TR<<"\n";
        ss<<"T1 = "; for(int i=0; i<n_T12; i++) ss<<T1[i]<<' '; ss<<'\n';
        ss<<"T2 = "; for(int i=0; i<n_T12; i++) ss<<T2[i]<<' '; ss<<'\n';
        ss<<"TE = "; for(int i=0; i<n_TE; i++) ss<<TE[i]*dt<<' '; ss<<'\n';

        ss<<"RF flip-angle   = "; for(int i=0; i<n_RF; i++) ss<<RF_FA[i]<<' '; ss<<'\n';
        ss<<"RF phase        = "; for(int i=0; i<n_RF; i++) ss<<RF_PH[i]<<' '; ss<<'\n';
        ss<<"RF time         = "; for(int i=0; i<n_RF; i++) ss<<RF_ST[i]*dt<<' '; ss<<'\n';

        ss<<"dephasing deg.  = "; for(int i=0; i<n_dephasing; i++) ss<<dephasing[i]<<' '; ss<<'\n';
        ss<<"dephasing time  = "; for(int i=0; i<n_dephasing; i++) ss<<dephasing_T[i]*dt<<' '; ss<<'\n';
        ss<<"gradient (x,y,z)=\n"; for(int i=0; i<n_gradient; i++) ss<<gradient_xyz[3*i+0]<<' '<<gradient_xyz[3*i+1]<<' '<<gradient_xyz[3*i+2]<<'\n';
        ss<<"gradient time   = "; for(int i=0; i<n_gradient; i++) ss<<gradient_T[i]*dt<<' '; ss<<'\n';
        ss<<"Cross Tissue Probability =\n"; for(int i=0; i<n_tissue_type; i++) {for(int j=0; j<n_tissue_type; j++) ss<<pXY[j+i*n_tissue_type]<<' '; ss<<'\n';};

        ss<<"sample length   = "<< sample_length[0] << " x " << sample_length[1] << " x " << sample_length[2] << " m" << '\n';
        ss<<"scale2grid      = "<< scale2grid[0] << " x " << scale2grid[1] << " x " << scale2grid[2] << '\n';
        ss<<"fieldmap size   = "<< fieldmap_size[0] << " x " << fieldmap_size[1] << " x " << fieldmap_size[2] << '\n';
        ss<<"matrix length   = "<< matrix_length << '\n';
        ss<<"diffusion const = "<<diffusion_const<<'\n'<<"dummy scans = "<<n_dummy_scan<<'\n'<<"spins = "<<n_spins<<'\n';
        ss<<"samples scales  = "<<n_sample_length_scales<<'\n'<<"timepoints = "<<n_timepoints<<'\n'<<"fieldmaps = "<<n_fieldmaps<<'\n';
        ss<<"max iterations  = "<<max_iterations<<'\n';
        ss<<"Tissue Types = " <<n_tissue_type << '\n'<< "Boundry Condition = " << enCrossBoundry << '\n';
        ss<<"Phase cycling   = "<<phase_cycling<<'\n'<<"Seed = "<<seed<<'\n' << "Record Trajectory = " << enRecordTrajectory << '\n';
        ss<<"Required CPU memory = "<<get_required_memory(1, "cpu")<<" MB\n";
        ss<<"Required GPU memory = "<<get_required_memory(1, "gpu")<<" MB\n";

        return ss.str();
    }

    size_t get_required_memory(uint8_t n_device=1, std::string type="gpu")
    {
        std::transform(type.begin(), type.end(), type.begin(), [](unsigned char c){ return std::tolower(c); }); 
        size_t spin = (type == "gpu") ? n_spins/n_device : n_spins;
        // fieldmap and mask
        size_t data_size_MB = fieldmap_size[0] * fieldmap_size[1] * fieldmap_size[2] * (sizeof(float) + sizeof(uint8_t)) / B2MB;
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

    void prepare()
    {
        c = cosf(RF_FA[0] * DEG2RAD); c2 = cosf(RF_FA[0] * DEG2RAD / 2.0f); 
        s = sinf(RF_FA[0] * DEG2RAD); s2 = sinf(RF_FA[0] * DEG2RAD / 2.0f);
        matrix_length = fieldmap_size[0] * fieldmap_size[1] * fieldmap_size[2];
        n_timepoints = TR / dt;
    }
} simulation_parameters;


inline void print_logo()
{ 
 std::cout << " \n"
" ____            _          __        __          _   _        \n"
"/ ___|   _ __   (_)  _ __   \\ \\      / /   __ _  | | | | __    \n"
"\\___ \\  | '_ \\  | | | '_ \\   \\ \\ /\\ / /   / _` | | | | |/ /    \n"
" ___) | | |_) | | | | | | |   \\ V  V /   | (_| | | | |   <     \n"
"|____/  | .__/  |_| |_| |_|    \\_/\\_/     \\__,_| |_| |_|\\_\\    \n"
"        |_|                                                    \n\n";

std::cout << "SpinWalk ver. " << SPINWALK_VERSION_MAJOR << "." << SPINWALK_VERSION_MINOR << "." << SPINWALK_VERSION_PATCH << std::endl;
}

#endif // __MISCELLANEOUS_H__