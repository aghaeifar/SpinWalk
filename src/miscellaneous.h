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

#include <cuda_runtime.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <vector>

#define DEG2RAD 0.0174532925199433 // = M_PI/180 
#define RAD2DEG 57.2957795130823


#define ERR_MSG  "\033[1;31mError:\033[0m "
#define ROUND(x) ((long)((x)+0.5))
#define MAX_RF 128          // maximum number of RF
#define MAX_TE 256          // maximum number of echo times
#define MAX_DEPHASE 128     // maximum number of dephasing

typedef struct simulation_parameters
{
    float T1, T2, TR, dt, B0, e1, e12, e2, e22, c, s, c2, s2;
    float RF_FA[MAX_RF], RF_PH[MAX_RF]; // refocusing FA
    float dephasing[MAX_DEPHASE]; // dephasing in degree
    int32_t RF_ST[MAX_RF], TE[MAX_TE], dephasing_T[MAX_DEPHASE]; // refocusing time in dt, echo times in dt, dephasing time in dt
    float sample_length[3], scale2grid[3], diffusion_const, phase_cycling;
    int32_t n_timepoints, n_sample_length_scales, n_fieldmaps, n_TE, n_RF, n_dephasing;
    int32_t n_dummy_scan;
    uint32_t n_spins, fieldmap_size[3], seed;
    uint64_t matrix_length;
    bool enDebug, enApplyFA2, enCrossBoundry;
    simulation_parameters():T1(2.2),T2(0.04),TR(0.04),dt(5e-5),B0(9.4),n_TE(0),n_RF(0),n_dummy_scan(0),phase_cycling(0.),enApplyFA2(false),enDebug(false),enCrossBoundry(true)
    {
        memset(fieldmap_size, 0, 3*sizeof(fieldmap_size[0])); 
        memset(sample_length, 0, 3*sizeof(sample_length[0]));
        memset(TE, 0, MAX_TE*sizeof(TE[0]));
        memset(RF_FA, 0, MAX_RF*sizeof(RF_FA[0]));
        memset(RF_ST, 0, MAX_RF*sizeof(RF_ST[0]));
        memset(RF_PH, 0, MAX_RF*sizeof(RF_PH[0]));
        memset(dephasing, 0, MAX_DEPHASE*sizeof(dephasing[0]));
        memset(dephasing_T, 0, MAX_DEPHASE*sizeof(dephasing_T[0]));
    }

    void dump()
    {
        std::cout<<"T1="<<T1<<" T2="<<T2<<" TR="<<TR<<" dt="<<dt<<" B0="<<B0<<'\n';
        std::cout<<"TE = "; for(int i=0; i<n_TE; i++) std::cout<<TE[i]*dt<<' '; std::cout<<'\n';
        std::cout<<"RF degree = "; for(int i=0; i<n_RF; i++) std::cout<<RF_FA[i]<<' '; std::cout<<'\n';
        std::cout<<"RF phase = "; for(int i=0; i<n_RF; i++) std::cout<<RF_PH[i]<<' '; std::cout<<'\n';
        std::cout<<"RF time = "; for(int i=0; i<n_RF; i++) std::cout<<RF_ST[i]*dt<<' '; std::cout<<'\n';
        std::cout<<"dephasing degree = "; for(int i=0; i<n_dephasing; i++) std::cout<<dephasing[i]<<' '; std::cout<<'\n';
        std::cout<<"dephasing time = "; for(int i=0; i<n_dephasing; i++) std::cout<<dephasing_T[i]*dt<<' '; std::cout<<'\n';
        std::cout<<"sample length = "<< sample_length[0] << " x " << sample_length[1] << " x " << sample_length[2] << " m" << '\n';
        std::cout<<"scale2grid = "<< scale2grid[0] << " x " << scale2grid[1] << " x " << scale2grid[2] << '\n';
        std::cout<<"fieldmap size = "<< fieldmap_size[0] << " x " << fieldmap_size[1] << " x " << fieldmap_size[2] << '\n';
        std::cout<<"diffusion const = "<<diffusion_const<<'\t'<<"dummy scans = "<<n_dummy_scan<<'\t'<<"spins = "<<n_spins<<'\n';
        std::cout<<"samples scales = "<<n_sample_length_scales<<'\t'<<"timepoints = "<<n_timepoints<<'\t'<<"fieldmaps = "<<n_fieldmaps<<'\n';
        std::cout<<"Apply FA/2 = "<<enApplyFA2<<'\t'<<"Boundry Condition = " << enCrossBoundry << '\n';
        std::cout<<"Phase cycling = "<<phase_cycling<<'\t'<<"Seed = "<<seed<<'\n';
        std::cout<<'\n';

        uint16_t fieldmap_size_MB = fieldmap_size[0] * fieldmap_size[1] * fieldmap_size[2] * (sizeof(float) + sizeof(char)) / 1024 / 1024;
        uint16_t variables_size_MB = n_spins * 3 *  (4 + n_TE) * sizeof(float) / 1024 / 1024;
        std::cout<<"Required GPU memory ≈ " << fieldmap_size_MB << " MB + " << variables_size_MB << " MB (fieldmap + variables)" << '\n';
        std::cout<<"Required RAM ≈ " << fieldmap_size_MB << " MB + " << variables_size_MB * n_sample_length_scales << " MB (fieldmap + variables)" << '\n';
    }

    void prepare()
    {
        c = cosf(RF_FA[0] * DEG2RAD); c2 = cosf(RF_FA[0] * DEG2RAD / 2.0f); 
        s = sinf(RF_FA[0] * DEG2RAD); s2 = sinf(RF_FA[0] * DEG2RAD / 2.0f);
        e1  = exp(-TR / T1); e12 = exp(-TR / (2. * T1));
        e2  = exp(-TR / T2); e22 = exp(-TR / (2. * T2));
        matrix_length = fieldmap_size[0] * fieldmap_size[1] * fieldmap_size[2];
        n_timepoints = TR / dt;
    }
} simulation_parameters;



#endif // __MISCELLANEOUS_H__