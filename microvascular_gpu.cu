/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: microvascular_gpu.cu
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

// compile :  nvcc microvascular_gpu.cu -Xptxas -v -O3  -arch=compute_86 -code=sm_86  -Xcompiler -fopenmp

#include <random>
#include <filesystem>

#include "./common/kernels.h"
#include "./common/miscellaneous.h"


#define CONFIG_FILE     "../inputs/config.ini"
#define THREADS_PER_BLOCK       64

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

    // ========== checking GPU(s) ==========
    int32_t device_count;
    checkCudaErrors(cudaGetDeviceCount(&device_count));
    std::cout << "Number of GPU(s): " << device_count << std::endl;
    param.n_spins /= device_count; // spins will be distributed in multiple GPUs (if there is). We hope it is divisible 

    // ========== load field-maps ==========
    pFieldMap = new float[param.matrix_length];
    std::vector<float> M0(3*param.n_spins*param.n_sample_length_scales*device_count, 0.f);
    std::vector<float> M1(3*param.n_spins*param.n_sample_length_scales*device_count, 0.f);

    for (int16_t fieldmap_no=0; fieldmap_no<param.n_fieldmaps; fieldmap_no++)
    {
        std::cout << "Loading fieldmap " << fieldmap_no+1 << " = " << filenames.at("fieldmap")[fieldmap_no] << std::endl;
        std::ifstream in_field(filenames.at("fieldmap")[fieldmap_no], std::ios::in | std::ios::binary);
        in_field.seekg(sizeof(int) * 3 + sizeof(float) * 3); // skip header for now, but should match with the mask
        in_field.read((char *)&pFieldMap[0], sizeof(float) * param.matrix_length);
        in_field.close();
    
        // ========== distributing between devices ==========
        cudaEvent_t start_0, end_0;
        checkCudaErrors(cudaEventCreate (&start_0));
        checkCudaErrors(cudaEventCreate (&end_0));
        checkCudaErrors(cudaEventRecord (start_0));
        #pragma omp parallel for
        for(int32_t d=0; d<device_count; d++)
        {
            cudaStream_t stream;
            checkCudaErrors(cudaSetDevice(d));            
            checkCudaErrors(cudaStreamCreate(&stream));

            simulation_parameters *d_param, param_local;
            float *d_pFieldMap, *d_position_start, *d_position_start_scaled;
            float*d_M1;
            bool *d_pMask;

            checkCudaErrors(cudaMalloc((void**)&d_param,                sizeof(simulation_parameters)));
            checkCudaErrors(cudaMalloc((void**)&d_pFieldMap,            sizeof(float) * param.matrix_length));    
            checkCudaErrors(cudaMalloc((void**)&d_position_start, sizeof(float) * 3 * param.n_spins));
            checkCudaErrors(cudaMalloc((void**)&d_position_start_scaled, sizeof(float) * 3 * param.n_spins));
            checkCudaErrors(cudaMalloc((void**)&d_M1,                 sizeof(float) * 3 * param.n_spins));
            checkCudaErrors(cudaMalloc((void**)&d_pMask,                sizeof(bool) * param.matrix_length));
            
            checkCudaErrors(cudaMemcpyAsync(d_pFieldMap, pFieldMap, sizeof(float) * param.matrix_length, cudaMemcpyHostToDevice, stream));
            checkCudaErrors(cudaMemcpyAsync(d_pMask, pMask, sizeof(bool) * param.matrix_length, cudaMemcpyHostToDevice, stream));
            checkCudaErrors(cudaMemcpyAsync(d_param, &param, sizeof(simulation_parameters), cudaMemcpyHostToDevice, stream));
            memcpy(&param_local, &param, sizeof(simulation_parameters));
            //std::cout<< "Memory allocated on GPU " << d << std::endl;

            // ========== run ==========
            int32_t numBlocks = (param.n_spins + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
            float elapsedTime;
            cudaEvent_t start;
            cudaEvent_t end;
            checkCudaErrors(cudaEventCreate (&start));
            checkCudaErrors(cudaEventCreate (&end));
            checkCudaErrors(cudaEventRecord (start));
            // generate initial spatial position for spins, based on sample_length_ref
            //std::cout<< d << ") Generating random initial position for spins... " ;
            generate_initial_position<<<numBlocks, THREADS_PER_BLOCK, 0, stream>>>(d_position_start, d_param, d_pMask);
            gpuCheckKernelExecutionError( __FILE__, __LINE__);
            //std::cout << "Done!" << std::endl;
            
            for (int32_t sl = 0; sl < param.n_sample_length_scales; sl++)
            {        
                // if (param.n_sample_length_scales > 1) 
                //     printf("%d, %2d) Simulating sample scale = %8.5f\n", d, sl, sample_length_scales[sl]);

                for(int i=0; i<3; i++)
                {
                    param_local.sample_length[i] = sample_length_scales[sl] * param.sample_length[i];    
                    param_local.scale2grid[i] = (param_local.fieldmap_size[i] - 1.) / param_local.sample_length[i];
                }   
                cudaMemcpy(d_param, &param_local, sizeof(simulation_parameters), cudaMemcpyHostToDevice);

                scale_initial_positions<<<numBlocks, THREADS_PER_BLOCK, 0, stream>>>(d_position_start_scaled, d_position_start, sample_length_scales[sl], param.n_spins);
                gpuCheckKernelExecutionError( __FILE__, __LINE__);

                simulation_kernel<<<numBlocks, THREADS_PER_BLOCK, 0, stream>>>(d_param, 
                                                                               d_pFieldMap, 
                                                                               d_pMask, 
                                                                               d_position_start_scaled, 
                                                                               d_M1);
                gpuCheckKernelExecutionError( __FILE__, __LINE__);                                                       
                checkCudaErrors(cudaMemcpyAsync(M1.data() + 3*param.n_spins*sl + 3*param.n_spins*param.n_sample_length_scales*d
                                                , d_M1, sizeof(float)*3*param.n_spins, cudaMemcpyDeviceToHost, stream));
            }    

            checkCudaErrors(cudaEventRecord(end));
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, end));
            printf("%d) Simulation took = %.2f seconds\n", d, elapsedTime/1000.);

            // ========== clean up GPU ==========
            checkCudaErrors(cudaFree(d_param));
            checkCudaErrors(cudaFree(d_pFieldMap));
            checkCudaErrors(cudaFree(d_pMask));
            checkCudaErrors(cudaFree(d_position_start));
            checkCudaErrors(cudaFree(d_position_start_scaled));
            checkCudaErrors(cudaFree(d_M1));
            checkCudaErrors(cudaEventDestroy(start));
            checkCudaErrors(cudaEventDestroy(end));
            checkCudaErrors(cudaStreamDestroy(stream));
        }

        float elapsedTime;
        checkCudaErrors(cudaEventRecord(end_0));
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start_0, end_0));
        std::cout << "Entire simulation over " << device_count << " GPU(s) took " << std::fixed << std::setprecision(2) << elapsedTime/1000. << " second(s)" << std::endl;
        checkCudaErrors(cudaEventDestroy(start_0));
        checkCudaErrors(cudaEventDestroy(end_0));
        
        // ========== save results ========== 
        std::string append = std::filesystem::path(filenames.at("fieldmap")[fieldmap_no]).stem().string(); // Thanks to C++17, we can use std::filesystem
        output_header hdr(3, param.n_spins, param.n_sample_length_scales, device_count);
        save_output(M1, filenames.at("output")[0], append, hdr, sample_length_scales);
        std::cout << std::string(50, '=') << std::endl;
    }

    // ========== clean up CPU ==========
    delete[] pFieldMap;
    delete[] pMask;
}

