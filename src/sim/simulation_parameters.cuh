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
#include <sstream>
#include <cmath>
#include <vector>
#include <boost/log/trivial.hpp> 

#ifdef __CUDACC__
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#endif

#define DEG2RAD 0.0174532925199433 // = M_PI/180 
#define RAD2DEG 57.2957795130823


// host vectors
typedef struct parameters_hvec
{
    std::vector<double> diffusivity;    
    std::vector<float> RF_FA_deg, RF_PH_deg, dephasing_deg, gradient_mTm, pXY, T1_ms, T2_ms;
    std::vector<int32_t> TE_us, RF_us, dephasing_us, gradient_us;
} parameters_hvec;

#ifdef __CUDACC__
// device vectors
typedef struct parameters_dvec
{
    thrust::device_vector<double> diffusivity;    
    thrust::device_vector<float> RF_FA_deg, RF_PH_deg, dephasing_deg, gradient_mTm, pXY, T1_ms, T2_ms;
    thrust::device_vector<int32_t> TE_us, RF_us, dephasing_us, gradient_us;
    void copy_from_host(const parameters_hvec &param_hvec)
    {
        diffusivity = param_hvec.diffusivity;
        RF_FA_deg = param_hvec.RF_FA_deg;
        RF_PH_deg = param_hvec.RF_PH_deg;
        dephasing_deg = param_hvec.dephasing_deg;
        gradient_mTm = param_hvec.gradient_mTm;
        pXY = param_hvec.pXY;
        T1_ms = param_hvec.T1_ms;
        T2_ms = param_hvec.T2_ms;
        TE_us = param_hvec.TE_us;
        RF_us = param_hvec.RF_us;
        dephasing_us = param_hvec.dephasing_us;
        gradient_us = param_hvec.gradient_us;
    }
} parameters_dvec;
#endif

template <typename T>
struct uvec
{
    const T *ptr = nullptr;
    size_t size = 0;
};

// universal vectors
typedef struct parameters_uvec
{
    uvec<double> diffusivity;    
    uvec<float> RF_FA_deg, RF_PH_deg, dephasing_deg, gradient_mTm, pXY, T1_ms, T2_ms;
    uvec<int32_t> TE_us, RF_us, dephasing_us, gradient_us;
    void copy_from_host(const parameters_hvec &param_hvec)
    {
        diffusivity.ptr = param_hvec.diffusivity.data(); 
        RF_FA_deg.ptr = param_hvec.RF_FA_deg.data();
        RF_PH_deg.ptr = param_hvec.RF_PH_deg.data();
        dephasing_deg.ptr = param_hvec.dephasing_deg.data();
        gradient_mTm.ptr = param_hvec.gradient_mTm.data();
        pXY.ptr = param_hvec.pXY.data();
        T1_ms.ptr = param_hvec.T1_ms.data();
        T2_ms.ptr = param_hvec.T2_ms.data();
        TE_us.ptr = param_hvec.TE_us.data();
        RF_us.ptr = param_hvec.RF_us.data();
        dephasing_us.ptr = param_hvec.dephasing_us.data();
        gradient_us.ptr = param_hvec.gradient_us.data();
        // fill the size
        diffusivity.size = param_hvec.diffusivity.size();
        RF_FA_deg.size = param_hvec.RF_FA_deg.size();
        RF_PH_deg.size = param_hvec.RF_PH_deg.size();
        dephasing_deg.size = param_hvec.dephasing_deg.size();
        gradient_mTm.size = param_hvec.gradient_mTm.size();
        pXY.size = param_hvec.pXY.size();
        T1_ms.size = param_hvec.T1_ms.size();
        T2_ms.size = param_hvec.T2_ms.size();
        TE_us.size = param_hvec.TE_us.size();
        RF_us.size = param_hvec.RF_us.size();
        dephasing_us.size = param_hvec.dephasing_us.size();
        gradient_us.size = param_hvec.gradient_us.size();
    }
#ifdef __CUDACC__
    void copy_from_device(const parameters_dvec &param_dvec)
    {
        diffusivity.ptr = thrust::raw_pointer_cast(param_dvec.diffusivity.data());
        RF_FA_deg.ptr = thrust::raw_pointer_cast(param_dvec.RF_FA_deg.data());
        RF_PH_deg.ptr = thrust::raw_pointer_cast(param_dvec.RF_PH_deg.data());
        dephasing_deg.ptr = thrust::raw_pointer_cast(param_dvec.dephasing_deg.data());
        gradient_mTm.ptr = thrust::raw_pointer_cast(param_dvec.gradient_mTm.data());
        pXY.ptr = thrust::raw_pointer_cast(param_dvec.pXY.data());
        T1_ms.ptr = thrust::raw_pointer_cast(param_dvec.T1_ms.data());
        T2_ms.ptr = thrust::raw_pointer_cast(param_dvec.T2_ms.data());
        TE_us.ptr = thrust::raw_pointer_cast(param_dvec.TE_us.data());
        RF_us.ptr = thrust::raw_pointer_cast(param_dvec.RF_us.data());
        dephasing_us.ptr = thrust::raw_pointer_cast(param_dvec.dephasing_us.data());
        gradient_us.ptr = thrust::raw_pointer_cast(param_dvec.gradient_us.data());
        // fill the size
        diffusivity.size = param_dvec.diffusivity.size();
        RF_FA_deg.size = param_dvec.RF_FA_deg.size();
        RF_PH_deg.size = param_dvec.RF_PH_deg.size();
        dephasing_deg.size = param_dvec.dephasing_deg.size();
        gradient_mTm.size = param_dvec.gradient_mTm.size();
        pXY.size = param_dvec.pXY.size();
        T1_ms.size = param_dvec.T1_ms.size();
        T2_ms.size = param_dvec.T2_ms.size();
        TE_us.size = param_dvec.TE_us.size();
        RF_us.size = param_dvec.RF_us.size();
        dephasing_us.size = param_dvec.dephasing_us.size();
        gradient_us.size = param_dvec.gradient_us.size();
    }
#endif
} parameters_uvec;


typedef struct parameters
{
    double fov[3];
    float B0=9.4, c, s;
    float phase_cycling;
    int32_t timestep_us, TR_us, n_dummy_scan=0;
    uint32_t n_spins=1e3, n_timepoints=0, n_substrate=0, n_scales=1;
    size_t phantom_size[3], seed=0, max_iterations=9999;
    int64_t matrix_length;
    bool enCrossFOV, enRecordTrajectory, fieldmap_exist;

    parameters() :
        TR_us(40e3),
        timestep_us(20),
        phase_cycling(0.),
        enCrossFOV(true),
        enRecordTrajectory(false),
        fieldmap_exist(true),
        matrix_length(0),
        n_spins(1e3)
    {
        memset(phantom_size,   0, 3*sizeof(phantom_size[0])); 
        memset(fov,   0, 3*sizeof(fov[0]));
    }
    

    // std::string dump()
    // {
    //     std::stringstream ss;
    //     ss<<"B0 = "<<B0<<" T\n";
    //     ss<<"timestep = "<<timestep_us<<" us.\n";
    //     ss<<"TR = "<<TR_us/1000.<<" ms.\n";
    //     ss<<"TE = "; for(int i=0; i<TE_us.data.size(); i++) ss<<TE_us.data[i]*timestep_us/1000.<<' '; ss<<"ms.\n";
    //     ss<<"T2 = "; for(int i=0; i<T2_ms.data.size(); i++) ss<<T2_ms.data[i]<<' '; ss<<"ms.\n";
    //     ss<<"T1 = "; for(int i=0; i<T1_ms.data.size(); i++) ss<<T1_ms.data[i]<<' '; ss<<"ms.\n";
    //     ss<<"sqrt(2*diffusivity*timestep) = "; for(int i=0; i<diffusivity.data.size(); i++) ss<<diffusivity.data[i]<<' '; ss<<"\n";
    //     ss<<"Cross Tissue Probability =\n"; for(int i=0; i<n_substrate; i++) {for(int j=0; j<n_substrate; j++) ss<<pXY.data[j+i*n_substrate]<<' '; ss<<'\n';};
    //     ss<<"RF flip-angle   = "; for(int i=0; i<RF_FA_deg.data.size(); i++) ss<<RF_FA_deg.data[i]<<' '; ss<<"deg.\n";
    //     ss<<"RF phase        = "; for(int i=0; i<RF_PH_deg.data.size(); i++) ss<<RF_PH_deg.data[i]<<' '; ss<<"deg.\n";
    //     ss<<"RF time         = "; for(int i=0; i<RF_us.data.size(); i++) ss<<RF_us.data[i]*timestep_us<<' '; ss<<"us.\n";
    //     ss<<"dephasing       = "; for(int i=0; i<dephasing_deg.data.size(); i++) ss<<dephasing_deg.data[i]<<' '; ss<<"deg.\n";
    //     ss<<"dephasing time  = "; for(int i=0; i<dephasing_us.data.size(); i++) ss<<dephasing_us.data[i]*timestep_us<<' '; ss<<"us.\n";
    //     ss<<"FoV             = "<< fov[0] << " x " << fov[1] << " x " << fov[2] << " m" << '\n';
    //     ss<<"fieldmap size   = "<< phantom_size[0] << " x " << phantom_size[1] << " x " << phantom_size[2] << '\n';
    //     ss<<"matrix length   = "<< matrix_length << '\n';
    //     ss<<"dummy scans     = "<< n_dummy_scan<<'\n';
    //     ss<<"spins           = "<< n_spins<<'\n';
    //     ss<<"timepoints      = "<< n_timepoints<<'\n';
    //     ss<<"max iterations  = "<< max_iterations<<'\n';
    //     ss<<"Pass FoV        = "<< enCrossFOV << '\n';
    //     ss<<"Phase cycling   = "<< phase_cycling<<'\n';
    //     ss<<"Seed            = "<< seed<<'\n';
    //     ss<<"off-resonance exists   = "<< (fieldmap_exist ? "Yes" : "No") <<'\n';
    //     ss<<"gradient (x,y,z) mT/m  =\n"; for(int i=0; i<gradient_mTm.data.size()/3; i++) ss<<gradient_mTm.data[3*i+0]<<' '<<gradient_mTm.data[3*i+1]<<' '<<gradient_mTm.data[3*i+2]<<'\n';
    //     ss<<"gradient  time         = "; for(int i=0; i<gradient_us.data.size(); i++) ss<<gradient_us.data[i]*timestep_us<<' '; ss<<"us\n";
    //     ss<<"Record Trajectory      = "<< (enRecordTrajectory ? "Yes" : "No")<<'\n';
    //     ss<<"Number of scales   = "<< n_scales<<'\n';
    //     ss<<"Number of substrates   = "<< n_substrate;
    //     return ss.str();
    // }

    bool prepare(parameters_hvec &param_hvec)
    {
        c = cosf(param_hvec.RF_FA_deg[0] * DEG2RAD); 
        s = sinf(param_hvec.RF_FA_deg[0] * DEG2RAD); 
        n_timepoints = TR_us / timestep_us;

        if (seed == 0)
            seed = std::random_device{}();

        for (int i = 0; i < n_substrate; i++)
            param_hvec.diffusivity[i] = 1e-3 * sqrt(2. * param_hvec.diffusivity[i] * timestep_us); // 1e-3 is the conversion factor from us to seconds for timestep_us

        if (n_dummy_scan < 0)
            n_dummy_scan = 5.0 * param_hvec.T1_ms[0] / TR_us * 1e3;

        return true;
    }
} parameters;

#endif // __SIMULATION_PARAMETERS_H__
