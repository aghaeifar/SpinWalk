
/* --------------------------------------------------------------------------
 * Project: 
 * File: 
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : 
 * -------------------------------------------------------------------------- */

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "ini.h"

#define ROUND(x) ((long)((x)+0.5))

#ifndef M_PI
#define M_PI 3.14159265359
#endif

typedef struct simulation_parameters
{
    float T1, T2, FA, TE, TR, dt, B0, e1, e12, e2, e22, c, s, c2, s2;
    float sample_length, diffusion_const, scale2grid;
    int16_t n_dummy_scan, n_timepoints, n_sample_length, n_fieldmaps;
    int32_t n_spins, fieldmap_size[3], seed, gpu_device_id;
    int64_t matrix_length;
    bool enDebug, enSteadyStateSimulation, enRefocusing180;
    void dump()
    {
        std::cout<<"T1="<<T1<<'\t'<<"T2="<<T2<<'\t'<<"FA="<<FA<<'\t'<<"TE="<<TE<<'\t'<<"TR="<<TR<<'\t'<<"dt="<<dt<<'\t'<<"B0="<<B0<<'\n';
        std::cout<<"sample length="<<sample_length<<'\t'<<"diffusion const="<<diffusion_const<<'\t'<<"dummy scans="<<n_dummy_scan<<'\t'<<"spins="<<n_spins<<'\n';
        std::cout<<"samples="<<n_sample_length<<'\t'<<"timepoints="<<n_timepoints<<'\t'<<"fieldmaps="<<n_fieldmaps<<'\t'<<"fieldmap size="<<fieldmap_size[0]<<" x "<<fieldmap_size[1]<<" x "<<fieldmap_size[2]<<'\n';
    }
#ifdef __CUDACC__
    __device__ void ddump()
    {
        printf("T1=%f n_spins=%d sample_length=%f sin=%f %f, cos=%f %f e1=%f %f e2=%f %f\n", T1, n_spins, sample_length, s, s2, c, c2, e1, e12, e2, e22);
    }
#endif
    void setFA(float flip_angle)
    {
        FA = flip_angle;
        c = cosf(FA); c2 = cosf(FA/2.); 
        s = sinf(FA); s2 = sinf(FA/2.);
    }
    void setTRT1T2(float tr, float t1=-1, float t2=-1)
    {
        TR = tr;
        T1 = t1>0 ? t1:T1;
        T2 = t2>0 ? t2:T2;
           
        e1  = exp(-TR / T1); e12 = exp(-TR / (2. * T1));
        e2  = exp(-TR / T2); e22 = exp(-TR / (2. * T2));
    }
} simulation_parameters;

typedef struct output_header
{
    int32_t n_spins, n_fieldmaps, n_sample_length, device_count;
    output_header(int32_t a, int32_t b, int32_t c, int32_t d): n_spins(a), n_fieldmaps(b), n_sample_length(c), device_count(d){}
} output_header;

typedef struct input_header
{
    int32_t size[3];
    float sample_length;
    input_header(int32_t *a, float b): sample_length(b){memcpy(size, a, 3*sizeof(int32_t));}
} input_header;

bool file_exist(std::string fname)
{
    std::ifstream f(fname.c_str());
    if(f.good() == false)
    {
        std::cout << "File does not exist: " << fname << std::endl;
        return false;
    }
    return true;
}

bool read_config(std::string config_filename, simulation_parameters &param, std::vector<float> &sample_length, std::map<std::string, std::vector<std::string> > &filenames)
{
    if(file_exist(config_filename.c_str()) == false)
        return false;
    // first, create a file instance
    mINI::INIFile file(config_filename);
    // next, create a structure that will hold data
    mINI::INIStructure ini;
    if (file.read(ini) == false)
    {
        std::cout << "Problem reading config file: " << config_filename << std::endl;
        return false;
    }
    if (ini.has("files") == false || ini.has("TISSUE_PARAMETERS") == false || ini.has("SIMULATION_PARAMETERS") == false || ini.has("SCAN_PARAMETERS") == false)
    {
        std::cout << "No all required sections in config file exist.";
        return false;
    }

    filenames.at("fieldmap").clear();
    for (int i = 0; ini.get("FILES").has("FIELD_MAP[" + std::to_string(i) + "]"); i++)
    {
        filenames.at("fieldmap").push_back(ini.get("FILES").get("FIELD_MAP[" + std::to_string(i) + "]"));
        if(file_exist(filenames.at("fieldmap").back().c_str()) == false)
            return false;
    }

    param.n_fieldmaps = filenames.at("fieldmap").size();

    filenames.at("mask").push_back(ini.get("FILES").get("MASK"));
    if(file_exist(filenames.at("mask").back().c_str()) == false)
        return false;
    filenames.at("output").push_back(ini.get("FILES").get("OUTPUTS"));

    param.dt = std::stof(ini.get("SCAN_PARAMETERS").get("DWELL_TIME"));
    param.n_dummy_scan = std::stoi(ini.get("SCAN_PARAMETERS").get("DUMMY_SCAN"));
    param.setTRT1T2(std::stof(ini.get("SCAN_PARAMETERS").get("TR")),
                    std::stof(ini.get("TISSUE_PARAMETERS").get("T1")),
                    std::stof(ini.get("TISSUE_PARAMETERS").get("T2"))
                    );
    param.TE = param.TR/2.;// std::stof(ini.get("SCAN_PARAMETERS").get("TE"));
    param.setFA(std::stof(ini.get("SCAN_PARAMETERS").get("FA")) * M_PI / 180.); // convert to radian

    param.B0 = std::stof(ini.get("SIMULATION_PARAMETERS").get("B0"));
    param.seed = std::stoi(ini.get("SIMULATION_PARAMETERS").get("SEED"));
    //en_phase_cycling = std::stoi(ini.get("SIMULATION_PARAMETERS").get("ENABLE_PHASE_CYCLING")) > 0;
    //en_180_refocusing = std::stoi(ini.get("SIMULATION_PARAMETERS").get("ENABLE_180_REFOCUSING")) > 0;
    param.n_spins = std::stof(ini.get("SIMULATION_PARAMETERS").get("NUMBER_OF_SPINS"));
    param.diffusion_const = std::stof(ini.get("SIMULATION_PARAMETERS").get("DIFFUSION_CONSTANT"));
    param.enRefocusing180 = ini.get("SIMULATION_PARAMETERS").get("ENABLE_180_REFOCUSING").compare("0") != 0;

    for (int i = 0; ini.get("SIMULATION_PARAMETERS").has("SAMPLE_LENGTH[" + std::to_string(i) + "]"); i++)
        sample_length.push_back(std::stof(ini.get("SIMULATION_PARAMETERS").get("SAMPLE_LENGTH[" + std::to_string(i) + "]")));
    param.n_sample_length = sample_length.size();

    param.gpu_device_id = std::stoi(ini.get("DEBUG").get("GPU_DEVICE_ID"));
    param.enDebug = ini.get("DEBUG").get("DUMP_INFO").compare("0") != 0;
    param.enSteadyStateSimulation = ini.get("DEBUG").get("SIMULATE_STEADYSTATE").compare("0") != 0;
    return true;
}


#ifdef __CUDACC__
__host__  __device__ __forceinline__
#else
inline
#endif
int32_t sub2ind(int32_t x, int32_t y, int32_t z, int32_t lenx, int32_t leny)
{
    return (z*lenx*leny + y*lenx + x);
}
