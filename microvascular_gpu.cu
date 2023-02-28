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
#include "./common/config_reader.h"


#define CONFIG_DEFAULT     "./inputs/config_default.ini"
#define THREADS_PER_BLOCK  64

using namespace std;

int main(int argc, char * argv[])
{
    std::vector<std::string> config_files(1, CONFIG_DEFAULT);
    if(argc > 2)
    {
        std::cout << "Usage: " << argv[0] << " <config_file>" << std::endl;
        return 1;
    }
    if(argc == 2)
        config_files.push_back(argv[1]);

    map<string, vector<string> > filenames = {{"fieldmap", vector<string>()},
                                              {"output", vector<string>()},
                                              {"XYZ0", vector<string>()},
                                              {"M0", vector<string>()} }; 
    std::vector<float> sample_length_scales, fieldmap;
    std::vector<char> mask;
    simulation_parameters param;
    input_header hdr_in;

    // ========== read config file ==========
    param.fieldmap_size[0] = param.fieldmap_size[1] = param.fieldmap_size[2] = 0;
    param.sample_length[0] = param.sample_length[1] = param.sample_length[2] = 0.f;
    for(uint8_t cnf_fl=0; cnf_fl<config_files.size(); cnf_fl++)
        if(read_config(config_files[cnf_fl], param, sample_length_scales, filenames) == false)
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

    // ========== Dump Settings ==========
    if(param.enDebug)
    {
        std::cout << "Dumping settings:" << std::endl;
        for (int32_t i = 0; i < param.n_fieldmaps; i++)
            std::cout << "Fieldmap " << i << " = " << filenames.at("fieldmap")[i] << std::endl;
        
        for (int32_t i = 0; i < param.n_sample_length_scales; i++)
            std::cout << "Sample length scale " << i << " = " << sample_length_scales[i] << std::endl;

        param.dump();
        std::cout<< std::string(30, '-')  << std::endl;
    }

    // ========== checking GPU(s) ==========
    int32_t device_count;
    checkCudaErrors(cudaGetDeviceCount(&device_count));
    std::cout << "Number of GPU(s): " << device_count << std::endl;
    param.n_spins /= device_count; // spins will be distributed in multiple GPUs (if there is). We hope it is divisible 

    std::vector<float> M0(3 * param.n_spins * param.n_sample_length_scales * device_count, 0.f);
    std::vector<float> M1(3 * param.n_spins * param.n_sample_length_scales * device_count, 0.f);
    std::vector<float> XYZ0(3 * param.n_spins * param.n_sample_length_scales * device_count, 0.f);
    std::vector<float> XYZ1(3 * param.n_spins * param.n_sample_length_scales * device_count, 0.f);

    for (int16_t fieldmap_no=0; fieldmap_no<param.n_fieldmaps; fieldmap_no++)
    {
        // ========== load field-maps ==========
        std::string fieldmap_file = filenames.at("fieldmap")[fieldmap_no];
        std::cout << "Loading fieldmap " << fieldmap_no << " = " << fieldmap_file << std::endl;
        std::ifstream in_field(fieldmap_file, std::ios::in | std::ios::binary);
        in_field.read((char *)&hdr_in, sizeof(input_header));
        std::copy(hdr_in.fieldmap_size, hdr_in.fieldmap_size+3, param.fieldmap_size);
        std::copy(hdr_in.sample_length, hdr_in.sample_length+3, param.sample_length);
        param.matrix_length = param.fieldmap_size[0] * param.fieldmap_size[1] * param.fieldmap_size[2];
        if(fieldmap.size() != param.matrix_length)
        {
            std::cout << "Fieldmap size changed. Re-allocating memory..." << std::endl;
            std::cout << "Old size: " << fieldmap.size() << std::endl;
            std::cout << "New size: " << param.matrix_length << std::endl;
            std::cout << "New length (um): " << param.sample_length[0]*1e6 << " " << param.sample_length[1]*1e6 << " " << param.sample_length[2]*1e6 << std::endl;
            fieldmap.resize(param.matrix_length);
            mask.resize(param.matrix_length);
        }
        in_field.read((char *)fieldmap.data(), sizeof(float) * param.matrix_length);
        in_field.read((char *)mask.data(),     sizeof(bool)  * param.matrix_length);
        in_field.close();

        for(int i=0; i<3; i++)
            param.scale2grid[i] = (param.fieldmap_size[i] - 1.) / param.sample_length[i];

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
            float *d_M1, *d_XYZ1;
            bool  *d_pMask;

            checkCudaErrors(cudaMalloc((void**)&d_param,                sizeof(simulation_parameters)));
            checkCudaErrors(cudaMalloc((void**)&d_pFieldMap,            fieldmap.size() * sizeof(fieldmap[0])));   
            checkCudaErrors(cudaMalloc((void**)&d_pMask,                mask.size() * sizeof(mask[0]))); 
            checkCudaErrors(cudaMalloc((void**)&d_position_start,       sizeof(float) * param.n_spins * 3));
            checkCudaErrors(cudaMalloc((void**)&d_position_start_scaled,sizeof(float) * param.n_spins * 3));
            checkCudaErrors(cudaMalloc((void**)&d_M1,                   sizeof(float) * param.n_spins * 3));
            checkCudaErrors(cudaMalloc((void**)&d_XYZ1,                 sizeof(float) * param.n_spins * 3));
            
            checkCudaErrors(cudaMemcpyAsync(d_pFieldMap, fieldmap.data(), fieldmap.size() * sizeof(fieldmap[0]), cudaMemcpyHostToDevice, stream));
            checkCudaErrors(cudaMemcpyAsync(d_pMask,     mask.data(),     mask.size() * sizeof(mask[0]),         cudaMemcpyHostToDevice, stream));
            checkCudaErrors(cudaMemcpyAsync(d_param,     &param,          sizeof(simulation_parameters),         cudaMemcpyHostToDevice, stream));
            memcpy(&param_local, &param, sizeof(simulation_parameters));

            // ========== run ==========
            int32_t numBlocks = (param.n_spins + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
            float elapsedTime;
            cudaEvent_t start;
            cudaEvent_t end;
            checkCudaErrors(cudaEventCreate (&start));
            checkCudaErrors(cudaEventCreate (&end));
            checkCudaErrors(cudaEventRecord (start));

            // generate initial spatial position for spins, based on sample_length_ref
            printf("GPU %d) Generating random initial position for spins... ", d);
            generate_initial_position<<<numBlocks, THREADS_PER_BLOCK, 0, stream>>>(d_position_start, d_param, d_pMask);
            gpuCheckKernelExecutionError( __FILE__, __LINE__);
            printf("Done!\n");
            
            for (int32_t sl = 0; sl < param.n_sample_length_scales; sl++)
            {        
                if (param.n_sample_length_scales > 1) 
                    printf("GPU %d) Simulating sample scale %2d = %8.5f\n", d, sl, sample_length_scales[sl]);

                for(int i=0; i<3; i++)
                {
                    param_local.sample_length[i] = sample_length_scales[sl] * param.sample_length[i];    
                    param_local.scale2grid[i] = (param_local.fieldmap_size[i] - 1.) / param_local.sample_length[i];
                }   
                cudaMemcpy(d_param, &param_local, sizeof(simulation_parameters), cudaMemcpyHostToDevice);

                scale_initial_positions<<<numBlocks, THREADS_PER_BLOCK, 0, stream>>>(d_position_start_scaled, d_position_start, sample_length_scales[sl], param.n_spins);
                gpuCheckKernelExecutionError( __FILE__, __LINE__);

                simulation_kernel<<<numBlocks, THREADS_PER_BLOCK, 0, stream>>>(d_param, d_pFieldMap, d_pMask, d_position_start_scaled, d_M1, d_XYZ1);
                gpuCheckKernelExecutionError( __FILE__, __LINE__);    

                int shift = 3*param.n_spins*sl + 3*param.n_spins*param.n_sample_length_scales*d ;                                        
                checkCudaErrors(cudaMemcpyAsync(M1.data()   + shift, d_M1,   sizeof(float)*3*param.n_spins, cudaMemcpyDeviceToHost, stream));
                checkCudaErrors(cudaMemcpyAsync(XYZ1.data() + shift, d_XYZ1, sizeof(float)*3*param.n_spins, cudaMemcpyDeviceToHost, stream));
            }    

            checkCudaErrors(cudaEventRecord(end));
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, end));
            printf("%d) Simulation took = %.2f second(s)\n", d, elapsedTime/1000.);

            // ========== clean up GPU ==========
            checkCudaErrors(cudaFree(d_param));
            checkCudaErrors(cudaFree(d_pFieldMap));
            checkCudaErrors(cudaFree(d_pMask));
            checkCudaErrors(cudaFree(d_position_start));
            checkCudaErrors(cudaFree(d_position_start_scaled));
            checkCudaErrors(cudaFree(d_M1));
            checkCudaErrors(cudaFree(d_XYZ1));
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
        std::string append = std::to_string(fieldmap_no) + "_" + std::filesystem::path(fieldmap_file).stem().string(); // Thanks to C++17, we can use std::filesystem
        output_header hdr(3, param.n_spins, param.n_sample_length_scales, device_count);
        save_output(M1  , filenames.at("output")[0], "M1_"   + append, hdr, sample_length_scales);
        save_output(XYZ1, filenames.at("output")[0], "XYZ1_" + append, hdr, sample_length_scales);

        std::cout << std::string(50, '=') << std::endl;
    }
}

