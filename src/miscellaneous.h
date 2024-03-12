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
#include <cmath>

#define SPINWALK_VERSION_MAJOR 1
#define SPINWALK_VERSION_MINOR 4
#define SPINWALK_VERSION_PATCH 5

#define DEG2RAD 0.0174532925199433 // = M_PI/180 
#define RAD2DEG 57.2957795130823

#define ERR_MSG  "\033[1;31mError:\033[0m "
#define ROUND(x) ((long)((x)+0.5))
#define MAX_RF 256          // maximum number of RF
#define MAX_TE 256          // maximum number of echo times
#define MAX_T12 256         // maximum number of relaxation times
#define MAX_DEPHASE 256     // maximum number of dephasing
#define MAX_GRADIENT 256    // maximum number of gradient

typedef struct simulation_parameters
{
    float TR, dt, B0, c, s, c2, s2;
    float T1[MAX_T12], T2[MAX_T12];
    float RF_FA[MAX_RF], RF_PH[MAX_RF]; // refocusing FA
    float dephasing[MAX_DEPHASE]; // dephasing in degree
    float gradient_xyz[3*MAX_GRADIENT]; // gradient in T/m
    int32_t RF_ST[MAX_RF], TE[MAX_TE], dephasing_T[MAX_DEPHASE], gradient_T[MAX_GRADIENT]; // refocusing time in dt, echo times in dt, dephasing time in dt
    float sample_length[3], scale2grid[3], diffusion_const, phase_cycling;
    int32_t n_timepoints, n_sample_length_scales, n_fieldmaps, n_TE, n_RF, n_dephasing, n_gradient, n_T12;
    int32_t n_dummy_scan;
    uint32_t n_spins, fieldmap_size[3], seed;
    uint64_t matrix_length;
    bool enDebug, enCrossBoundry, enMultiTissue, enRecordTrajectory;
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
        phase_cycling(0.),
        enDebug(false),
        enCrossBoundry(true),
        enMultiTissue(false),
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

        std::fill(T1, T1 + MAX_T12, 2.2);
        std::fill(T2, T2 + MAX_T12, 0.04);
    }

    std::string dump()
    {
        std::stringstream ss;
        ss<<"TR="<<TR<<" dt="<<dt<<" B0="<<B0<<'\n';
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

        ss<<"sample length   = "<< sample_length[0] << " x " << sample_length[1] << " x " << sample_length[2] << " m" << '\n';
        ss<<"scale2grid      = "<< scale2grid[0] << " x " << scale2grid[1] << " x " << scale2grid[2] << '\n';
        ss<<"fieldmap size   = "<< fieldmap_size[0] << " x " << fieldmap_size[1] << " x " << fieldmap_size[2] << '\n';
        ss<<"diffusion const = "<<diffusion_const<<'\t'<<"dummy scans = "<<n_dummy_scan<<'\t'<<"spins = "<<n_spins<<'\n';
        ss<<"samples scales  = "<<n_sample_length_scales<<'\t'<<"timepoints = "<<n_timepoints<<'\t'<<"fieldmaps = "<<n_fieldmaps<<'\n';
        ss<<"Multi-Tissues   = "<<enMultiTissue<<'\t'<<"Boundry Condition = " << enCrossBoundry << '\n';
        ss<<"Phase cycling   = "<<phase_cycling<<'\t'<<"Seed = "<<seed<<'\n';
        ss<<'\n';

        auto fieldmap_size_MB = fieldmap_size[0] * fieldmap_size[1] * fieldmap_size[2] * (sizeof(float) + sizeof(char)) / 1024 / 1024;
        auto variables_size_MB = n_spins * 3 *  (4 + n_TE) * sizeof(float) / 1024 / 1024;
        ss<<"Required GPU memory ≈ " << fieldmap_size_MB << " MB + " << variables_size_MB << " MB (fieldmap + variables)" << '\n';
        ss<<"Required RAM ≈ " << fieldmap_size_MB << " MB + " << variables_size_MB * n_sample_length_scales << " MB (fieldmap + variables)" << '\n';

        return ss.str();
    }

    uint32_t get_required_memory(uint32_t n_device=1)
    {
        auto fieldmap_mask_size_MB = fieldmap_size[0] * fieldmap_size[1] * fieldmap_size[2] * (sizeof(float) + sizeof(char)) / 1024 / 1024;
        auto variables_size_MB = n_spins * 3 *  (4 + n_TE) * sizeof(float) / 1024 / 1024 / n_device; // 3 for x,y,z, 4 for M0, XYZ0, XYZ0_scaled, XYZ1,  n_TE for M1
        auto trajectory_size_MB = enRecordTrajectory ? n_spins * n_timepoints * 3 * sizeof(float) / 1024 / 1024 / n_device: 0; // 3 for x,y,z,
        return fieldmap_mask_size_MB + variables_size_MB + trajectory_size_MB;
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