/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: spinwalk.cu
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating randomwalk in microvascular network
 * -------------------------------------------------------------------------- */

// compile(lin) :  nvcc ./src/spinwalk.cu ./src/kernels.cu -I ./include/ -Xptxas -v -O3  -arch=compute_75 -code=sm_75  -Xcompiler -fopenmp -o spinwalk
// compile(win) :  nvcc ./src/spinwalk.cu ./src/kernels.cu -I ./include/ -Xptxas -v -O3  -arch=compute_86 -code=sm_86  -Xcompiler /openmp -std=c++17 -o spinwalk

#include <chrono>
#include <random>
#include <iomanip>
#include <filesystem>
#include <cuda_runtime.h>
#include "tqdm.h"
#include "kernels.cuh"
#include "file_utils.h"
#include "helper_cuda.h"
#include "miscellaneous.h"
#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

#define THREADS_PER_BLOCK  64
#define LOG_FILE "spinwalk.log"

namespace bl = boost::log;

bool simulate(simulation_parameters param, std::map<std::string, std::vector<std::string> > filenames, std::vector<double> sample_length_scales)
{
    // ========== checking number of GPU(s) ==========
    int32_t device_count;
    checkCudaErrors(cudaGetDeviceCount(&device_count));
    BOOST_LOG_TRIVIAL(info) << "Number of available GPU(s): " << device_count; 
    
    param.n_spins /= device_count; // spins will be distributed in multiple GPUs (if there is). We hope it is divisible 
    size_t numBlocks = (param.n_spins + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    
    // ========== allocate memory on CPU ==========
    auto st = std::chrono::steady_clock::now();
    size_t trj  = param.enRecordTrajectory ? param.n_timepoints * (param.n_dummy_scan + 1) : 1;
    size_t len0 = 3 * param.n_spins * device_count;
    size_t len1 = len0 * param.n_sample_length_scales * trj;
    size_t len2 = len0 * param.n_sample_length_scales * param.n_TE;
    BOOST_LOG_TRIVIAL(info) << "Memory size: XYZ0=" << len0 << ", XYZ1=" << len1 << ", M0=" << len0 << ", M1=" << len2 << "";
    std::vector<float> fieldmap(param.matrix_length, 0.f);
    std::vector<uint8_t> mask(param.matrix_length, 0);
    std::vector<float> XYZ0(len0, 0.f);     // memory layout(column-wise): [3 x n_spins]
    std::vector<float> XYZ1(len1, 0.f);     // memory layout(column-wise): [3 x timepoints x n_spins x n_sample_length_scales] or [3 x 1 x n_spins x n_sample_length_scales]
    std::vector<float> M0(len0, 0.f);       // memory layout(column-wise): [3 x n_spins]
    std::vector<float> M1(len2, 0.f);       // memory layout(column-wise): [3 x n_TE x n_spins x n_sample_length_scales]
    std::vector<uint8_t> T(M1.size()/3, 0); // memory layout(column-wise): [1 x n_TE x n_spins x n_sample_length_scales]
    BOOST_LOG_TRIVIAL(info) << "Memory allocation (CPU) took " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - st).count() << " ms";
    // ========== allocate memory on GPU ==========
    std::vector<float *> d_pFieldMap(device_count, NULL);
    std::vector<float *> d_M0(device_count, NULL), d_M1(device_count, NULL);
    std::vector<float *> d_XYZ1(device_count, NULL), d_XYZ0(device_count, NULL), d_XYZ0_scaled(device_count, NULL);
    std::vector<uint8_t *>  d_T(device_count, NULL);
    std::vector<uint8_t *>  d_pMask(device_count, NULL);
    std::vector<simulation_parameters *> d_param(device_count, NULL);
    std::vector<cudaStream_t> streams(device_count, NULL);
    
    #pragma omp parallel for
    for(int32_t d=0; d<device_count; d++)
    {  
        checkCudaErrors(cudaSetDevice(d));            
        checkCudaErrors(cudaStreamCreate(&streams[d]));
        // allocate memory on GPU
        checkCudaErrors(cudaMalloc((void**)&d_param[d],         sizeof(simulation_parameters)));
        checkCudaErrors(cudaMalloc((void**)&d_pFieldMap[d],     sizeof(fieldmap[0]) * fieldmap.size()));   
        checkCudaErrors(cudaMalloc((void**)&d_pMask[d],         sizeof(mask[0]) * mask.size())); 
        checkCudaErrors(cudaMalloc((void**)&d_XYZ0[d],          sizeof(XYZ0[0]) * XYZ0.size() / device_count));
        checkCudaErrors(cudaMalloc((void**)&d_XYZ0_scaled[d],   sizeof(XYZ0[0]) * XYZ0.size() / device_count));
        checkCudaErrors(cudaMalloc((void**)&d_XYZ1[d],          sizeof(XYZ1[0]) * XYZ1.size() / device_count / param.n_sample_length_scales));
        checkCudaErrors(cudaMalloc((void**)&d_M0[d],            sizeof(M0[0]) * M0.size() / device_count));
        checkCudaErrors(cudaMalloc((void**)&d_M1[d],            sizeof(M1[0]) * M1.size() / device_count / param.n_sample_length_scales));
        checkCudaErrors(cudaMalloc((void**)&d_T[d],             sizeof(T[0]) * T.size() / device_count / param.n_sample_length_scales));
        checkCudaErrors(cudaStreamSynchronize(streams[d]));    
    }
    cudaEvent_t start;
    cudaEvent_t end;  
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&end));

    std::cout << std::string(50, '=') << std::endl;
    for (int16_t fieldmap_no=0; fieldmap_no<param.n_fieldmaps; fieldmap_no++)
    {
        bool hasXYZ0 = false;
        // ========== load files (field-maps, xyz0, m0) ==========
        if(file_utils::read_fieldmap(filenames.at("FIELDMAP")[fieldmap_no], fieldmap, mask, param) == false)
            return false;

        if(fieldmap_no < filenames.at("XYZ0").size())
        {   
            if(file_utils::read_file(filenames.at("XYZ0")[fieldmap_no], XYZ0) == false)
                return false;
            hasXYZ0 = true;
        }

        if(fieldmap_no < filenames.at("M0").size())
        {
            if(file_utils::read_file(filenames.at("M0")[fieldmap_no], M0) == false)
                return false;
        }
        else
        {   // all spins are aligned with B0 (M0 = (0, 0, 1))
            long index = 0;
            BOOST_LOG_TRIVIAL(info) << "Generating M0(0, 0, 1)..." << std::endl;
            std::generate(M0.begin(), M0.end(), [&index](){return (index++ % 3 == 2) ? 1.f : 0.f;});
        }

        for(int i=0; i<M0.size()/3; i += M0.size()/3/2)
            BOOST_LOG_TRIVIAL(info) << "M0 of the spin " << i << " = (" << M0[3*i] << ", " << M0[3*i+1] << ", " << M0[3*i+2] << ")" << std::endl;
        
        for(int i=0; i<3; i++)
            param.scale2grid[i] = param.fieldmap_size[i] / param.sample_length[i];
        
        if (hasXYZ0 && param.n_sample_length_scales > 1)
        {
            BOOST_LOG_TRIVIAL(error) << "loading XYZ0 from file while having more than 1 sample length scales is not supported!" << std::endl;
            return false;
        }
        
        // ========== copy data to GPU(s) ==========       
        #pragma omp parallel for
        for(int32_t d=0; d<device_count; d++)
        {  
            simulation_parameters param_local;
            memcpy(&param_local, &param, sizeof(simulation_parameters));
            param_local.seed += d * param.n_spins; // different seed for each GPU
            checkCudaErrors(cudaSetDevice(d));            
            checkCudaErrors(cudaMemcpyAsync(d_pFieldMap[d], fieldmap.data(),        fieldmap.size()*sizeof(fieldmap[0]), cudaMemcpyHostToDevice, streams[d]));
            checkCudaErrors(cudaMemcpyAsync(d_pMask[d],     mask.data(),            mask.size() * sizeof(mask[0]),       cudaMemcpyHostToDevice, streams[d]));
            checkCudaErrors(cudaMemcpyAsync(d_param[d],     &param_local,           sizeof(simulation_parameters),       cudaMemcpyHostToDevice, streams[d]));
            checkCudaErrors(cudaMemcpyAsync(d_M0[d],        &M0[3*param.n_spins*d], 3*param.n_spins*sizeof(M0[0]),       cudaMemcpyHostToDevice, streams[d]));

            // convert fieldmap from Tesla to degree per dwell time
            cu_scaleArray<<<uint64_t(fieldmap.size()/THREADS_PER_BLOCK)+1, THREADS_PER_BLOCK, 0, streams[d]>>>(d_pFieldMap[d], param.B0*GAMMA*param.dt*RAD2DEG, fieldmap.size());
            if(hasXYZ0 == false)
            {   // generate initial spatial position for spins, based on sample_length_ref
                BOOST_LOG_TRIVIAL(info) << "GPU " << d << ") Generating random initial position for spins... (seed = " << param_local.seed << " grid = " << numBlocks << " x " << THREADS_PER_BLOCK << ")";
                cu_randPosGen<<<numBlocks, THREADS_PER_BLOCK, 0, streams[d]>>>(d_XYZ0[d], d_param[d], d_pMask[d]);
                gpuCheckKernelExecutionError( __FILE__, __LINE__);
                BOOST_LOG_TRIVIAL(info) << "GPU " << d << ") Done!";
            }
            else // copy initial spatial position and magnetization for spins
                checkCudaErrors(cudaMemcpyAsync(d_XYZ0[d], &XYZ0[3*param.n_spins*d], 3*param.n_spins*sizeof(XYZ0[0]), cudaMemcpyHostToDevice, streams[d]));  
            checkCudaErrors(cudaStreamSynchronize(streams[d]));    
        }

        // ========== run ==========       
        checkCudaErrors(cudaEventRecord(start));
        
        tqdm bar;
        simulation_parameters param_local;
        memcpy(&param_local, &param, sizeof(simulation_parameters));
        for (int32_t sl = 0; sl < param.n_sample_length_scales; sl++)
        {
            for (int i = 0; i < 3; i++)
            {
                param_local.sample_length[i] = sample_length_scales[sl] * param.sample_length[i];
                param_local.scale2grid[i] = param_local.fieldmap_size[i] / param_local.sample_length[i];
            }
            
            #pragma omp parallel for
            for (int32_t d = 0; d < device_count; d++)
            {   
                BOOST_LOG_TRIVIAL(info) << "GPU " << d << ") Fieldmap " << fieldmap_no << ", simulating sample length scale " << sample_length_scales[sl];     
                checkCudaErrors(cudaSetDevice(d));
                cudaMemcpy(d_param[d], &param_local, sizeof(simulation_parameters), cudaMemcpyHostToDevice);

                cu_scalePos<<<numBlocks, THREADS_PER_BLOCK, 0, streams[d]>>>(d_XYZ0_scaled[d], d_XYZ0[d], sample_length_scales[sl], param.n_spins);
                gpuCheckKernelExecutionError(__FILE__, __LINE__);
                
                cu_sim<<<numBlocks, THREADS_PER_BLOCK, 0, streams[d]>>>(d_param[d], d_pFieldMap[d], d_pMask[d], d_M0[d], d_XYZ0_scaled[d], d_M1[d], d_XYZ1[d], d_T[d]); // d_XYZ0_scaled[d]
                gpuCheckKernelExecutionError(__FILE__, __LINE__);

                size_t shift = 3*param.n_TE*param.n_spins*device_count*sl + 3*param.n_TE*param.n_spins*d;
                checkCudaErrors(cudaMemcpyAsync(M1.data()   + shift, d_M1[d]  , 3*param.n_TE*param.n_spins*sizeof(M1[0]), cudaMemcpyDeviceToHost, streams[d]));                
                shift = 3*param.n_spins*trj*device_count*sl + 3*param.n_spins*trj*d;
                checkCudaErrors(cudaMemcpyAsync(XYZ1.data() + shift, d_XYZ1[d], 3*param.n_spins*trj*sizeof(XYZ1[0]), cudaMemcpyDeviceToHost, streams[d]));
                shift = param.n_TE*param.n_spins*device_count*sl + param.n_TE*param.n_spins*d;
                checkCudaErrors(cudaMemcpyAsync(T.data() + shift, d_T[d], param.n_TE*param.n_spins*sizeof(T[0]), cudaMemcpyDeviceToHost, streams[d]));
            } 
            bar.progress(sl, param.n_sample_length_scales);
        }
        bar.finish();

        float elapsedTime;
        checkCudaErrors(cudaEventRecord(end));
        checkCudaErrors(cudaDeviceSynchronize());        
        checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, end));
        std::cout << "Simulation over " << device_count << " GPU(s) took " << std::fixed << std::setprecision(2) << elapsedTime/1000. << " second(s)" << '\n';

        // ========== save results ========== 
        std::cout << "Saving the results to disk" << '\n';
        file_utils::output_header hdr(3, param.n_TE, param.n_spins * device_count, param.n_sample_length_scales);
        file_utils::save_output((char*)M1.data(), M1.size()*sizeof(M1[0]),filenames.at("M1")[fieldmap_no], hdr, sample_length_scales);

        hdr.dim2 = param.enRecordTrajectory ? param.n_timepoints * (param.n_dummy_scan + 1) : 1;
        file_utils::save_output((char*)XYZ1.data(), XYZ1.size()*sizeof(XYZ1[0]), filenames.at("XYZ1")[fieldmap_no], hdr, sample_length_scales);

        hdr.dim1 = 1;
        hdr.dim2 = param.n_TE;
        file_utils::save_output((char*)T.data(), T.size()*sizeof(T[0]), filenames.at("T")[fieldmap_no], hdr, sample_length_scales);

        std::cout << std::string(50, '=') << std::endl;
    }

    // ========== clean up GPU ==========
    #pragma omp parallel for
    for(int32_t d=0; d<device_count; d++)
    {
        checkCudaErrors(cudaSetDevice(d));   
        checkCudaErrors(cudaFree(d_param[d]));
        checkCudaErrors(cudaFree(d_pFieldMap[d]));
        checkCudaErrors(cudaFree(d_pMask[d]));
        checkCudaErrors(cudaFree(d_XYZ0[d]));
        checkCudaErrors(cudaFree(d_XYZ0_scaled[d]));
        checkCudaErrors(cudaFree(d_M1[d]));
        checkCudaErrors(cudaFree(d_XYZ1[d]));
        checkCudaErrors(cudaStreamDestroy(streams[d]));            
    }
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(end));
    
    return true;
} 


bool dump_settings(simulation_parameters param, std::map<std::string, std::vector<std::string> > filenames, std::vector<double> sample_length_scales)
{
    std::stringstream ss;
    ss  << "Dumping settings:" << '\n';
    for (std::map<std::string, std::vector<std::string>>::iterator it=filenames.begin(); it!=filenames.end(); ++it)
        for (int i = 0; i< it->second.size(); i++)
            ss << it->first << "[" << i << "] = " << it->second.at(i) << '\n';
    
    ss<< "\nSample length scale = [";
    for (int32_t i = 0; i < param.n_sample_length_scales; i++)
        ss << sample_length_scales[i] << ", ";
    ss << "]\n";
    
    file_utils::input_header hdr_in;
    if(file_utils::read_header(filenames.at("FIELDMAP")[0], hdr_in) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "reading header of fieldmap " << filenames.at("FIELDMAP")[0] << " failed. Aborting...!";
        return false;
    }
    std::copy(hdr_in.fieldmap_size, hdr_in.fieldmap_size+3, param.fieldmap_size);
    std::copy(hdr_in.sample_length, hdr_in.sample_length+3, param.sample_length);
    ss << param.dump();
    BOOST_LOG_TRIVIAL(info) << ss.str();
    return true;
}

int main(int argc, char * argv[])
{
    bool bStatus = true;
    print_logo();
    // ========== parse command line arguments ==========
    boost::program_options::options_description desc("Options");
    desc.add_options()
        ("help,h", "help message (this menu)")
        ("sim_off,s", "no simulation, only read config files")
        ("configs,c", boost::program_options::value<std::vector<std::string>>()->multitoken(), "config. files as many as you want. e.g. -c config1.ini config2.ini ... configN.ini");

    boost::program_options::variables_map vm;
    boost::program_options::parsed_options parsed = boost::program_options::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
    boost::program_options::store(parsed, vm);    
    std::vector<std::string> unreg = boost::program_options::collect_unrecognized(parsed.options, boost::program_options::include_positional);
    boost::program_options::notify(vm);

    bl::add_file_log(bl::keywords::file_name=LOG_FILE, bl::keywords::target_file_name = LOG_FILE, bl::keywords::format = "[%TimeStamp%] [%Severity%]: %Message%", bl::keywords::auto_flush = true);
    bl::add_common_attributes();
    
    // ========== print help ==========
    if (vm.count("help") || vm.count("configs") == 0 || argc == 1 || unreg.size() > 0)
    {
        std::cout << desc;
        print_device_info();
        return 1;
    }

    std::vector<std::string> config_files = vm["configs"].as<std::vector<std::string>>();  
    bool bNoSim = vm.count("sim_off") > 0;

    // ========== loop over configs and simulate ==========
    std::cout << "Running simulation for " << config_files.size() << " config(s)..." << std::endl;
    auto start = std::chrono::steady_clock::now();
    for(const auto& cfile : config_files)
    {
        std::map<std::string, std::vector<std::string> > filenames = {{"FIELDMAP", 	std::vector<std::string>()},  // input:  map of off-resonance in Tesla
                                                                      {"XYZ0", 		std::vector<std::string>()},  // input:  spins starting spatial positions in meters
                                                                      {"XYZ1", 		std::vector<std::string>()},  // output: spins last spatial positions in meters
                                                                      {"M0", 		std::vector<std::string>()},  // input:  spins initial magnetization
                                                                      {"M1", 		std::vector<std::string>()},  // output: spins final magnetization
                                                                      {"T", 		std::vector<std::string>()}}; // output: tissue index

        std::vector<double> sample_length_scales;
        simulation_parameters param;
        
        // ========== read config file ==========
        bStatus &= file_utils::read_config(cfile, param, sample_length_scales, filenames);
        
        if (param.seed == 0)
            param.seed = std::random_device{}();

        param.n_timepoints = param.TR / param.dt; // includes start point
        
        // ========== dump settings ==========
        bStatus &= dump_settings(param, filenames, sample_length_scales);

        // ========== Check GPU memory ==========
        bStatus &= check_memory_size(param.get_required_memory(getDeviceCount()));

        if (bStatus == false)
        {
            std::cout << ERR_MSG << "Simulation failed. See the log file " << LOG_FILE <<", Aborting...!" << std::endl;
            return 1;
        }
        
        if (bNoSim)
            continue;
        std::cout << "Simulation starts..." << std::endl;
        if(simulate(param, filenames, sample_length_scales) == false)
        {
            std::cout << ERR_MSG << "Simulation failed. See the log file " << LOG_FILE <<", Aborting...!" << std::endl;
            return 1;
        }
    }
    std::cout << "Simulation(s) finished successfully! Entire elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()/1000. << " second(s)."  << std::endl;
    return 0;
}
