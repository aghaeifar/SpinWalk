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

#include <sstream>
#include <random>
#include <filesystem>
#include <execution>
#include <numeric> // std::inner_product, std::iota
#include <chrono>
#include <string.h>
#include "../common/kernels.h"
#include "../common/miscellaneous.h"

#define CONFIG_FILE     "../inputs/config.ini"

using namespace std;

int main(int argc, char * argv[])
{
    std::string config_file = CONFIG_FILE;
    if(argc > 2)
    {
        std::cout << "Usage: " << argv[0] << " <config_file>" << std::endl;
        return 1;
    }
    if(argc == 2)
        config_file = argv[1];

    map<string, vector<string> > filenames = {{"fieldmap", vector<string>()},
                                              {"mask", vector<string>()},
                                              {"output", vector<string>()} }; 
    std::vector<float> sample_length_scales;
    simulation_parameters param;
    float *pFieldMap = NULL;
    bool *pMask = NULL;

    // ========== read config file ==========
    if(read_config(config_file, param, sample_length_scales, filenames) == false)
    {
        std::cout << "Reading config file failed. Aborting...!" << std::endl;
        return 1;
    }

    if (param.seed == 0)
        param.seed = std::random_device{}();

    param.n_timepoints = param.TR / param.dt; // includes start point

    // ========== simulating steady-state signal ==========
    if(param.enSteadyStateSimulation && param.n_dummy_scan != 0)
    {
        simulate_steady_state(param);
        std::cout<< std::string(30, '-')  << std::endl;
    }

    // ========== load mask ==========
    // read mask
    std::ifstream in_mask(filenames.at("mask")[0], std::ios::in | std::ios::binary);
    in_mask.read((char *)&param.fieldmap_size[0], sizeof(int) * 3);
    in_mask.read((char *)&param.sample_length[0], sizeof(float) * 3);
    param.matrix_length = param.fieldmap_size[0] * param.fieldmap_size[1] * param.fieldmap_size[2];
    pMask = new bool[param.matrix_length];
    in_mask.read((char *)&pMask[0], sizeof(bool) * param.matrix_length);
    in_mask.close();

    for(int i=0; i<3; i++)
        param.scale2grid[i] = (param.fieldmap_size[i] - 1.) / param.sample_length[i];  
    
    // ========== Dump Settings ==========
    if(param.enDebug)
    {
        std::cout << "Dumping settings:" << std::endl;
        for (int32_t i = 0; i < param.n_fieldmaps; i++)
            std::cout << "Fieldmap " << i+1 << " = " << filenames.at("fieldmap")[i] << std::endl;
        
        for (int32_t i = 0; i < param.n_sample_length_scales; i++)
            std::cout << "Sample length scale " << i+1 << " = " << sample_length_scales[i] << std::endl;

        param.dump();
        std::cout<< std::string(30, '-')  << std::endl;
    }

    // ========== load field-maps ==========
    pFieldMap = new float[param.matrix_length];
    simulation_parameters param_local;
    std::vector<float> position_start(3*param.n_spins, 0.f);
    std::vector<float> position_end(position_start);
    std::vector<float> position_start_scaled(position_start);
    std::vector<float> M0(3*param.n_spins*param.n_sample_length_scales, 0.f);
    std::vector<float> M1(3*param.n_spins*param.n_sample_length_scales, 0.f);
    memcpy(&param_local, &param, sizeof(simulation_parameters));

    std::vector<int> a(param.n_spins);
    std::iota(a.begin(), a.end(), 0);

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

        auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start);
        std::cout << "Simulation took " << elapsedTime.count() << " second(s)" << std::endl;
        
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