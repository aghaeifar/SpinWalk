/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: microvascular_cpu.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

// compile :  g++ microvascular_cpu.cpp -ltbb -O3 -g -Wall -ftree-vectorizer-verbose=5 -msse -msse2 -msse3 -march=native -mtune=native -ffast-math

#undef __CUDACC__

#include <random> 
#include <execution>
#include <numeric> // std::inner_product, std::iota
#include <chrono>
#include <random>
#include <filesystem>

#include "./common/kernels.h"
#include "./common/reader.h"

#define CONFIG_FILE     "../inputs/config.ini"

using namespace std;

bool simulate(simulation_parameters param, std::map<std::string, std::vector<std::string> > filenames, std::vector<float> sample_length_scales)
{
    std::vector<float> fieldmap;
    std::vector<char> mask;
    simulation_parameters param_local;

    uint32_t len0 = 3 * param.n_spins;
    uint32_t len1 = len0 * param.n_sample_length_scales;
    std::vector<float> M0(len0, 0.f); 
    std::vector<float> M1(len1, 0.f);
    std::vector<float> XYZ0(len0, 0.f);
    std::vector<float> XYZ1(len1, 0.f);

    std::vector<int> a(param.n_spins);
    std::iota(a.begin(), a.end(), 0);

    std::cout << std::string(50, '=') << std::endl;
    for (int16_t fieldmap_no=0; fieldmap_no<param.n_fieldmaps; fieldmap_no++)
    {
        bool hasXYZ0 = false;
        // ========== load files (field-maps, xyz0, m0) ==========
        if(reader::read_fieldmap(filenames.at("fieldmap")[fieldmap_no], fieldmap, mask, param) == false)
            return false;

        if(filenames.at("xyz0")[fieldmap_no].empty() == false)
        {
            if(reader::read_file(filenames.at("xyz0")[fieldmap_no], XYZ0) == false)
                return false;
            
            std::cout << "Checking XYZ0 is not in the mask..." << std::endl;
            uint32_t t = is_masked(XYZ0, mask, &param);
            if(t>0)
            {
                std::cout << ERR_MSG << t << " element(s) of XYZ0 is in the mask. Aborting...!" << std::endl;
                return 1;
            }
            hasXYZ0 = true;
        }

        if(filenames.at("m0")[fieldmap_no].empty() == false)
        {
            if(reader::read_file(filenames.at("m0")[fieldmap_no], M0) == false)
                return false;
        }
        else
        {   // all spins are aligned with B0 (M0 = (0, 0, 1))
            long index = 0;
            std::cout << "Generating M0(0, 0, 1)..." << std::endl;
            std::generate(M0.begin(), M0.end(), [&index](){return (index++ % 3 == 2) ? 1.f : 0.f;});
        }

        if(param.enDebug)
            for(int i=0; i<M0.size()/3; i += M0.size()/3/2)
                std::cout << "M0 of the spin " << i << " = (" << M0[3*i] << ", " << M0[3*i+1] << ", " << M0[3*i+2] << ")" << std::endl;

        for(int i=0; i<3; i++)
            param.scale2grid[i] = (param.fieldmap_size[i] - 1.) / param.sample_length[i];
        
        if (hasXYZ0 && param.n_sample_length_scales > 1)
        {
            std::cout << ERR_MSG << "loading XYZ0 from file while having more than 1 sample length scales is not supported!" << std::endl;
            return false;
        }

        if(hasXYZ0 == false)
        {   // generate initial spatial position for spins, based on sample_length_ref
            std::for_each(std::execution::par_unseq, std::begin(a), std::end(a), [&](int spin_no)
            {
                generate_initial_position(XYZ0.data(), &param, (bool *)mask.data(), spin_no);
            });
        }

    // generate initial spatial position for spins, based on sample_length_ref
    std::cout << "Generating random initial position for spins... ";
    std::for_each(std::execution::par_unseq, std::begin(a), std::end(a), [&](int spin_no)
    {
        generate_initial_position(position_start.data(), &param, pMask, spin_no);
    });

    std::cout << "Done!" << std::endl;

    for (int16_t fieldmap_no=0; fieldmap_no<param.n_fieldmaps; fieldmap_no++)
    {
        std::cout << "Loading fieldmap " << fieldmap_no+1 << " = " << filenames.at("fieldmap")[fieldmap_no] << std::endl;
        std::ifstream in_field(filenames.at("fieldmap")[fieldmap_no], std::ios::in | std::ios::binary);
        in_field.seekg(sizeof(int) * 3 + sizeof(float) * 3); // skip header for now, but should match with the mask
        in_field.read((char *)&pFieldMap[0], sizeof(float) * param.matrix_length);
        in_field.close();       

        // ========== run ==========
        auto start = std::chrono::system_clock::now();
        for (int32_t sl = 0; sl < param.n_sample_length_scales; sl++)
        {
            if (param.n_sample_length_scales > 1)
                printf("%2d / %2d) Simulating sample scale = %8.5f\n", sl, param.n_sample_length_scales, sample_length_scales[sl]);

            for (int i = 0; i < 3; i++)
            {
                param_local.sample_length[i] = sample_length_scales[sl] * param.sample_length[i];
                param_local.scale2grid[i] = (param_local.fieldmap_size[i] - 1.) / param_local.sample_length[i];
            }
            std::for_each(std::execution::par_unseq, std::begin(a), std::end(a), [&](int spin_no)
            {
                for(int i=0; i<3; i++)
                    position_start_scaled[3*spin_no + i] = position_start[3*spin_no + i] * sample_length_scales[sl];

                simulation_kernel(&param_local, pFieldMap, pMask, position_start_scaled.data(), M1.data() + 3*param.n_spins*sl, spin_no);
            });
        }    

        auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
        std::cout << "Simulation took " << std::fixed << std::setprecision(2) << elapsedTime.count()/1000. << " second(s)" << std::endl;
        
        // ========== save results ==========      
        std::string append = std::filesystem::path(filenames.at("fieldmap")[fieldmap_no]).stem().string(); // Thanks to C++17, we can use std::filesystem
        output_header hdr(3, param.n_spins, param.n_sample_length_scales);
        save_output(M1, filenames.at("output")[0], append, hdr, sample_length_scales);
        std::cout << std::string(50, '=') << std::endl;
    }

    // ========== clean up CPU ==========
    delete[] pFieldMap;
    delete[] pMask;
}

int main(int argc, char * argv[])
{
    std::vector<std::string> config_files;
    if(argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <config_file>" << std::endl;
        return 1;
    }
    for(uint8_t i=1; i<argc; i++)
        config_files.push_back(argv[i]);

    std::cout << "Running " << config_files.size() << " simulation(s)..." << std::endl;
    for(uint8_t cnf=0; cnf<config_files.size(); cnf++)
    {
        map<string, vector<string> > filenames = {{"fieldmap", vector<string>()},
                                                  {"xyz0", vector<string>()},
                                                  {"xyz1", vector<string>()},
                                                  {"m0", vector<string>()},
                                                  {"m1", vector<string>()} }; 

        std::vector<float> sample_length_scales;
        simulation_parameters param;

        // ========== read config file ==========
        param.fieldmap_size[0] = param.fieldmap_size[1] = param.fieldmap_size[2] = 0;
        param.sample_length[0] = param.sample_length[1] = param.sample_length[2] = 0.f;
        if(reader::read_config(config_files[cnf], param, sample_length_scales, filenames) == false)
        {
            std::cout << ERR_MSG << "reading config file failed. Aborting...!" << std::endl;
            return 1;
        }

        if (param.seed == 0)
            param.seed = std::random_device{}();

        param.n_timepoints = param.TR / param.dt; // includes start point

        // ========== simulating steady-state signal ==========
        if(param.enSteadyStateSimulation)
        {
            simulate_steady_state(param);
            std::cout<< std::string(30, '-')  << std::endl;
        }

        // ========== Dump Settings ==========
        if(param.enDebug)
        {
            std::cout << "Dumping settings:" << std::endl;
            for (std::map<std::string, std::vector<std::string>>::iterator it=filenames.begin(); it!=filenames.end(); ++it, std::cout << std::endl)
                for (int i = 0; i< it->second.size(); i++)
                    std::cout << it->first << "[" << i << "] = " << it->second.at(i) << std::endl;
            
            for (int32_t i = 0; i < param.n_sample_length_scales; i++)
                std::cout << "Sample length scale " << i << " = " << sample_length_scales[i] << std::endl;

            param.dump();
            std::cout<< std::string(30, '-')  << std::endl;
        }

        if(simulate(param, filenames, sample_length_scales) == false)
            return 1;
    }
    std::cout << "Simulation(s) finished successfully!" << std::endl;
    return 0;
}
