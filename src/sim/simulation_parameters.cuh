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
#include <boost/log/trivial.hpp> 
#include "host_device_array.cuh"

#define DEG2RAD 0.0174532925199433 // = M_PI/180 
#define RAD2DEG 57.2957795130823

typedef struct simulation_parameters
{
    double fov[3];
    float B0=9.4, c, s;
    float phase_cycling;
    int32_t timestep_us, TR_us, n_dummy_scan=0;
    uint32_t n_spins=1e3, n_timepoints=0, n_substrate=0, n_fov_scale=1;
    size_t phantom_size[3], seed=0, max_iterations=9999;
    int64_t matrix_length;
    bool enCrossFOV, enRecordTrajectory, fieldmap_exist;

    host_device_array<double> diffusivity;    
    host_device_array<float> RF_FA_deg, RF_PH_deg, dephasing_deg, gradient_mTm, pXY, T1_ms, T2_ms;
    host_device_array<int32_t> TE_us, RF_us, dephasing_us, gradient_us;
    
    simulation_parameters() :
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
    
    void init(bool device, bool copy_from_host = false)
    {
        diffusivity.init(device, copy_from_host);
        RF_FA_deg.init(device, copy_from_host);
        RF_PH_deg.init(device, copy_from_host);
        dephasing_deg.init(device, copy_from_host);
        gradient_mTm.init(device, copy_from_host);
        pXY.init(device, copy_from_host);
        T1_ms.init(device, copy_from_host);
        T2_ms.init(device, copy_from_host);
        TE_us.init(device, copy_from_host);
        RF_us.init(device, copy_from_host);
        dephasing_us.init(device, copy_from_host);
        gradient_us.init(device, copy_from_host);
    }

    std::string dump()
    {
        std::stringstream ss;
        ss<<"B0 = "<<B0<<" T\n";
        ss<<"timestep = "<<timestep_us<<" us.\n";
        ss<<"TR = "<<TR_us/1000.<<" ms.\n";
        ss<<"TE = "; for(int i=0; i<TE_us.data.size(); i++) ss<<TE_us.data[i]*timestep_us/1000.<<' '; ss<<"ms.\n";
        ss<<"T2 = "; for(int i=0; i<T2_ms.data.size(); i++) ss<<T2_ms.data[i]<<' '; ss<<"ms.\n";
        ss<<"T1 = "; for(int i=0; i<T1_ms.data.size(); i++) ss<<T1_ms.data[i]<<' '; ss<<"ms.\n";
        ss<<"sqrt(2*diffusivity*timestep) = "; for(int i=0; i<diffusivity.data.size(); i++) ss<<diffusivity.data[i]<<' '; ss<<"\n";
        ss<<"Cross Tissue Probability =\n"; for(int i=0; i<n_substrate; i++) {for(int j=0; j<n_substrate; j++) ss<<pXY.data[j+i*n_substrate]<<' '; ss<<'\n';};
        ss<<"RF flip-angle   = "; for(int i=0; i<RF_FA_deg.data.size(); i++) ss<<RF_FA_deg.data[i]<<' '; ss<<"deg.\n";
        ss<<"RF phase        = "; for(int i=0; i<RF_PH_deg.data.size(); i++) ss<<RF_PH_deg.data[i]<<' '; ss<<"deg.\n";
        ss<<"RF time         = "; for(int i=0; i<RF_us.data.size(); i++) ss<<RF_us.data[i]*timestep_us<<' '; ss<<"us.\n";
        ss<<"dephasing       = "; for(int i=0; i<dephasing_deg.data.size(); i++) ss<<dephasing_deg.data[i]<<' '; ss<<"deg.\n";
        ss<<"dephasing time  = "; for(int i=0; i<dephasing_us.data.size(); i++) ss<<dephasing_us.data[i]*timestep_us<<' '; ss<<"us.\n";
        ss<<"FoV             = "<< fov[0] << " x " << fov[1] << " x " << fov[2] << " m" << '\n';
        ss<<"fieldmap size   = "<< phantom_size[0] << " x " << phantom_size[1] << " x " << phantom_size[2] << '\n';
        ss<<"matrix length   = "<< matrix_length << '\n';
        ss<<"dummy scans     = "<< n_dummy_scan<<'\n';
        ss<<"spins           = "<< n_spins<<'\n';
        ss<<"timepoints      = "<< n_timepoints<<'\n';
        ss<<"max iterations  = "<< max_iterations<<'\n';
        ss<<"Pass FoV        = "<< enCrossFOV << '\n';
        ss<<"Phase cycling   = "<< phase_cycling<<'\n';
        ss<<"Seed            = "<< seed<<'\n';
        ss<<"off-resonance exists   = "<< (fieldmap_exist ? "Yes" : "No") <<'\n';
        ss<<"gradient (x,y,z) mT/m  =\n"; for(int i=0; i<gradient_mTm.data.size(); i++) ss<<gradient_mTm.data[3*i+0]<<' '<<gradient_mTm.data[3*i+1]<<' '<<gradient_mTm.data[3*i+2]<<'\n';
        ss<<"gradient  time         = "; for(int i=0; i<gradient_us.data.size(); i++) ss<<gradient_us.data[i]*timestep_us<<' '; ss<<"us\n";
        ss<<"Record Trajectory      = "<< (enRecordTrajectory ? "Yes" : "No")<<'\n';
        ss<<"Number of FoV scales   = "<< n_fov_scale<<'\n';
        ss<<"Number of substrates   = "<< n_substrate;
        return ss.str();
    }

    bool prepare()
    {
        c = cosf(RF_FA_deg.data[0] * DEG2RAD); 
        s = sinf(RF_FA_deg.data[0] * DEG2RAD); 
        n_timepoints = TR_us / timestep_us;

        if (seed == 0)
            seed = std::random_device{}();

        for (int i = 0; i < n_substrate; i++)
            diffusivity.data[i] = 1e-3 * sqrt(2. * diffusivity.data[i] * timestep_us); // 1e-3 is the conversion factor from us to seconds for timestep_us

        if (n_dummy_scan < 0)
            n_dummy_scan = 5.0 * T1_ms.data[0] / TR_us * 1e3;

        return true;
    }
} simulation_parameters;




#endif // __SIMULATION_PARAMETERS_H__