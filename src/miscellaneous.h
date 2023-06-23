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
#include <filesystem>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "rotation.h"


#define ERR_MSG  "\033[1;31mError:\033[0m "
#define ROUND(x) ((long)((x)+0.5))
#define MAX_SE 20   // maximum number of spin-echoes
#define MAX_TE 60   // maximum number of echo times


typedef struct simulation_parameters
{
    float T1, T2, FA, TR, dt, B0, e1, e12, e2, e22, c, s, c2, s2;
    float RF_SE[MAX_SE], RF_SE_PHS[MAX_SE]; // refocusing FA
    uint16_t T_SE[MAX_SE], TE[MAX_TE]; // refocusing time in dt, echo times in dt
    float sample_length[3], scale2grid[3], diffusion_const, phase_cycling;
    uint16_t n_dummy_scan, n_timepoints, n_sample_length_scales, n_fieldmaps, n_TE, n_SE;
    uint32_t n_spins, fieldmap_size[3], seed;
    uint64_t matrix_length;
    bool enDebug, enSteadyStateSimulation, enRefocusing, enApplyFA2;
    simulation_parameters():T1(2.2),T2(0.04),FA(16),TR(0.04),dt(5e-5),B0(9.4),n_TE(0),n_SE(0),n_dummy_scan(0),phase_cycling(0.),enSteadyStateSimulation(false),enRefocusing(false),enApplyFA2(false),enDebug(false)
    {
        memset(fieldmap_size, 0, 3*sizeof(fieldmap_size[0])); 
        memset(sample_length, 0, 3*sizeof(sample_length[0]));
        memset(TE, 0, MAX_TE*sizeof(TE[0]));
        memset(RF_SE, 0, MAX_SE*sizeof(RF_SE[0]));
        memset(T_SE, 0, MAX_SE*sizeof(T_SE[0]));
    }

    void dump()
    {
        std::cout<<"T1="<<T1<<" T2="<<T2<<" FA="<<FA<<" TR="<<TR<<" dt="<<dt<<" B0="<<B0<<'\n';
        std::cout<<"TE = "; for(int i=0; i<n_TE; i++) std::cout<<TE[i]*dt<<' '; std::cout<<'\n';
        std::cout<<"Refocusing RF degree = "; for(int i=0; i<n_SE; i++) std::cout<<RF_SE[i]<<' '; std::cout<<'\n';
        std::cout<<"Refocusing RF time = "; for(int i=0; i<n_SE; i++) std::cout<<T_SE[i]*dt<<' '; std::cout<<'\n';
        std::cout<<"sample length = "<< sample_length[0] << " x " << sample_length[1] << " x " << sample_length[2] << " m" << '\n';
        std::cout<<"scale2grid = "<< scale2grid[0] << " x " << scale2grid[1] << " x " << scale2grid[2] << '\n';
        std::cout<<"fieldmap size = "<< fieldmap_size[0] << " x " << fieldmap_size[1] << " x " << fieldmap_size[2] << '\n';
        std::cout<<"diffusion const = "<<diffusion_const<<'\t'<<"dummy scans = "<<n_dummy_scan<<'\t'<<"spins = "<<n_spins<<'\n';
        std::cout<<"samples scales = "<<n_sample_length_scales<<'\t'<<"timepoints = "<<n_timepoints<<'\t'<<"fieldmaps = "<<n_fieldmaps<<'\n';
        std::cout<<"Refocusing = "<<enRefocusing<<'\t'<<"Apply FA/2 = "<<enApplyFA2<<'\t'<<"Simulate steady-state = "<<enSteadyStateSimulation<<'\n';
        std::cout<<"Phase cycling = "<<phase_cycling<<'\t'<<"Seed = "<<seed<<'\n';
        std::cout<<'\n';

        uint16_t fieldmap_size_MB = fieldmap_size[0] * fieldmap_size[1] * fieldmap_size[2] * (sizeof(float) + sizeof(char)) / 1024 / 1024;
        uint16_t variables_size_MB = n_spins * 3 *  (4 + n_TE) * sizeof(float) / 1024 / 1024;
        std::cout<<"Required GPU memory ≈ " << fieldmap_size_MB << " MB + " << variables_size_MB << " MB (fieldmap + variables)" << '\n';
        std::cout<<"Required RAM ≈ " << fieldmap_size_MB << " MB + " << variables_size_MB * n_sample_length_scales << " MB (fieldmap + variables)" << '\n';

        size_t free, total;
        cudaMemGetInfo(&free, &total);
        std::cout << "Free GPU memory: " << free / 1024 / 1024 << " MB (out of " << total / 1024 / 1024 << " MB)" << std::endl;
    }
#ifdef __CUDACC__
    __device__ void ddump()
    {
        printf("T1=%f n_spins=%d sample_length=%f x %f x  %f sin=%f %f, cos=%f %f e1=%f %f e2=%f %f\n", T1, n_spins, sample_length[0], sample_length[1], sample_length[2], s, s2, c, c2, e1, e12, e2, e22);
    }
#endif

    void prepare()
    {
        c = cosf(FA * DEG2RAD); c2 = cosf(FA * DEG2RAD / 2.0f); 
        s = sinf(FA * DEG2RAD); s2 = sinf(FA * DEG2RAD / 2.0f);
        e1  = exp(-TR / T1); e12 = exp(-TR / (2. * T1));
        e2  = exp(-TR / T2); e22 = exp(-TR / (2. * T2));
        matrix_length = fieldmap_size[0] * fieldmap_size[1] * fieldmap_size[2];
        n_timepoints = TR / dt;
    }
} simulation_parameters;


typedef struct output_header
{
    int32_t dim1, dim2, dim3, dim4;
    output_header(int32_t a, int32_t b=1, int32_t c=1, int32_t d=1): dim1(a), dim2(b), dim3(c), dim4(d){}
} output_header;

typedef struct input_header
{
    uint32_t fieldmap_size[3];
    float sample_length[3];
    input_header(uint32_t *a, float *b) {memcpy(fieldmap_size, a, 3*sizeof(uint32_t)); memcpy(sample_length, b, 3*sizeof(float));}
    input_header(){};
    void print()
    {
        std::cout << "Size = " << fieldmap_size[0] << " x " << fieldmap_size[1] << " x " << fieldmap_size[2] << std::endl;
        std::cout << "Length = " << sample_length[0]*1e6 << " x " << sample_length[1]*1e6 << " x " << sample_length[2]*1e6 << " um^3" << std::endl;
    }
} input_header;


bool save_output(std::vector<float> &data, std::string output_filename, output_header hdr, std::vector<float> &additional_hdr)
{
    std::cout << "Saving output to: " << std::filesystem::absolute(output_filename) << std::endl;
    std::filesystem::path parent_path = std::filesystem::absolute(output_filename).parent_path();
    if (std::filesystem::is_directory(parent_path) == false)
    {
        std::cout << ERR_MSG << "cannot find directory " << parent_path.string() << ". Trying to create it." << std::endl;
        if(std::filesystem::create_directories(parent_path) == false)
        {
            std::cout << ERR_MSG << "cannot create directory " << parent_path.string() << std::endl;
            return false;
        }
    }

    std::ofstream file(output_filename, std::ios::out | std::ios::binary);
    if (file.is_open() == false)
    {
        std::cout << ERR_MSG << "cannot open file " << std::filesystem::absolute(output_filename) << std::endl;
        return false;
    }
    int32_t header_size = sizeof(output_header) + additional_hdr.size() * sizeof(additional_hdr[0]);
    file.write((char*)&header_size, sizeof(int32_t));
    file.write((char*)&hdr, sizeof(output_header));
    file.write((char*)additional_hdr.data(), additional_hdr.size() * sizeof(additional_hdr[0]));
    file.write((char*)data.data(), data.size() * sizeof(data[0]));
    file.close();
    return true; 
}

#endif // __MISCELLANEOUS_H__