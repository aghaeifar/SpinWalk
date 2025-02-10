/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: simulation_parameters.cuh
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : 
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


template <typename T>
std::string vec2str(const std::vector<T>& vec) {
    std::ostringstream oss;
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(oss, ", "));
    std::string content = oss.str();
    if (content.size() >= 2) 
        content.erase(content.size() - 2);
    return content;
}

// host vectors
typedef struct parameters_hvec
{
    std::vector<double> diffusivity;    
    std::vector<float> RF_FA_deg, RF_PH_deg, dephasing_deg, gradientX_mTm, gradientY_mTm, gradientZ_mTm, pXY, T1_ms, T2_ms;
    std::vector<int32_t> TE_us, RF_us, dephasing_us, gradient_us;

    std::string dump() 
    {
        std::stringstream ss;
        ss<<"TE = " << vec2str(TE_us) <<" us.\n";
        ss<<"T2 = " << vec2str(T2_ms) <<" ms.\n";
        ss<<"T1 = " << vec2str(T1_ms) <<" ms.\n";        
        ss<<"RF flip-angle   = " << vec2str(RF_FA_deg) <<" deg.\n";
        ss<<"RF phase        = " << vec2str(RF_PH_deg) <<" deg.\n";
        ss<<"RF time         = " << vec2str(RF_us) << " us.\n";
        ss<<"dephasing       = " << vec2str(dephasing_deg) << " deg.\n";
        ss<<"dephasing time  = " << vec2str(dephasing_us) << " us.\n";
        ss<<"gradient X      = " << vec2str(gradientX_mTm) << " mT/m\n";
        ss<<"gradient Y      = " << vec2str(gradientY_mTm) << " mT/m\n";
        ss<<"gradient Z      = " << vec2str(gradientZ_mTm) <<" mT/m\n";
        ss<<"gradient time   = " << vec2str(gradient_us) << " us\n";
        ss<<"diffusivity     = "<< vec2str(diffusivity) <<"\n";
        int n_substrate = sqrt(pXY.size());
        ss<<"Cross Tissue Probability =\n"; for(int i=0; i<n_substrate; i++) {for(int j=0; j<n_substrate; j++) ss<<pXY[j+i*n_substrate]<<' '; ss<<'\n';};
        return ss.str();
    }
} parameters_hvec;

#ifdef __CUDACC__
// device vectors
typedef struct parameters_dvec
{
    thrust::device_vector<double> diffusivity;    
    thrust::device_vector<float> RF_FA_deg, RF_PH_deg, dephasing_deg, gradientX_mTm, gradientY_mTm, gradientZ_mTm, pXY, T1_ms, T2_ms;
    thrust::device_vector<int32_t> TE_us, RF_us, dephasing_us, gradient_us;
    void copy_from_host(const parameters_hvec &param_hvec)
    {
        diffusivity = param_hvec.diffusivity;
        RF_FA_deg = param_hvec.RF_FA_deg;
        RF_PH_deg = param_hvec.RF_PH_deg;
        dephasing_deg = param_hvec.dephasing_deg;
        gradientX_mTm = param_hvec.gradientX_mTm;
        gradientY_mTm = param_hvec.gradientY_mTm;
        gradientZ_mTm = param_hvec.gradientZ_mTm;
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
    uvec<float> RF_FA_deg, RF_PH_deg, dephasing_deg, gradientX_mTm, gradientY_mTm, gradientZ_mTm, pXY, T1_ms, T2_ms;
    uvec<int32_t> TE_us, RF_us, dephasing_us, gradient_us;
    void copy_from_host(const parameters_hvec &param_hvec)
    {
        diffusivity.ptr = param_hvec.diffusivity.data(); 
        RF_FA_deg.ptr = param_hvec.RF_FA_deg.data();
        RF_PH_deg.ptr = param_hvec.RF_PH_deg.data();
        dephasing_deg.ptr = param_hvec.dephasing_deg.data();
        gradientX_mTm.ptr = param_hvec.gradientX_mTm.data();
        gradientY_mTm.ptr = param_hvec.gradientY_mTm.data();
        gradientZ_mTm.ptr = param_hvec.gradientZ_mTm.data();
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
        gradientX_mTm.size = param_hvec.gradientX_mTm.size();
        gradientY_mTm.size = param_hvec.gradientY_mTm.size();
        gradientZ_mTm.size = param_hvec.gradientZ_mTm.size();
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
        gradientX_mTm.ptr = thrust::raw_pointer_cast(param_dvec.gradientX_mTm.data());
        gradientY_mTm.ptr = thrust::raw_pointer_cast(param_dvec.gradientY_mTm.data());
        gradientZ_mTm.ptr = thrust::raw_pointer_cast(param_dvec.gradientZ_mTm.data());
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
        gradientX_mTm.size = param_dvec.gradientX_mTm.size();
        gradientY_mTm.size = param_dvec.gradientY_mTm.size();
        gradientZ_mTm.size = param_dvec.gradientZ_mTm.size();
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
        TR_us(-1),
        timestep_us(-1),
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
    

    std::string dump()
    {
        std::stringstream ss;
        ss<<"B0              = "<< B0<<" T\n";
        ss<<"timestep        = "<< timestep_us <<" us.\n";
        ss<<"TR              = "<< TR_us/1000. <<" ms.\n";
        ss<<"FoV             = "<< fov[0] << " x " << fov[1] << " x " << fov[2] << " m" << '\n';
        ss<<"fieldmap size   = "<< phantom_size[0] << " x " << phantom_size[1] << " x " << phantom_size[2] << '\n';
        ss<<"matrix length   = "<< matrix_length << '\n';
        ss<<"dummy scans     = "<< n_dummy_scan <<'\n';
        ss<<"spins           = "<< n_spins <<'\n';
        ss<<"timepoints      = "<< n_timepoints <<'\n';
        ss<<"max iterations  = "<< max_iterations <<'\n';
        ss<<"Pass FoV        = "<< enCrossFOV << '\n';
        ss<<"Phase cycling   = "<< phase_cycling <<'\n';
        ss<<"Seed            = "<< seed <<'\n';
        ss<<"off-resonance   = "<< (fieldmap_exist ? "Yes" : "No") <<'\n';
        ss<<"Save Trajectory = "<< (enRecordTrajectory ? "Yes" : "No")<<'\n';
        ss<<"N scales        = "<< n_scales <<'\n';
        ss<<"N substrates    = "<< n_substrate <<'\n';;
        return ss.str();
    }

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
