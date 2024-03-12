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

#include <random>
#include <filesystem>
#include <iomanip>
#include "helper_cuda.h"
#include <cuda_runtime.h>
#include "kernels.cuh"
#include "file_utils.h"
#include "miscellaneous.h"
#include "tqdm.h"

#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

#define THREADS_PER_BLOCK  64
#define LOG_FILE "spinwalk.log"

using namespace std;

bool simulate(simulation_parameters param, std::map<std::string, std::vector<std::string> > filenames, std::vector<float> sample_length_scales)
{
    std::vector<float> fieldmap;
    std::vector<uint8_t> mask;
    // ========== checking GPU(s) ==========
    int32_t device_count;
    checkCudaErrors(cudaGetDeviceCount(&device_count));

    param.n_spins /= device_count; // spins will be distributed in multiple GPUs (if there is). We hope it is divisible 
    size_t numBlocks = (param.n_spins + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

    size_t len0 = 3 * param.n_spins * device_count;
    size_t len1 = len0 * param.n_sample_length_scales;
    size_t len2 = len1 * param.n_TE;
    size_t len3 = len0 * param.n_timepoints * param.enRecordTrajectory;        
    std::vector<float> XYZ0(len0, 0.f); // memory layout(column-wise): [3 x n_spins]
    std::vector<float> XYZ1(len1, 0.f); // memory layout(column-wise): [3 x n_spins x n_sample_length_scales]
    std::vector<float> M0(len0, 0.f);   // memory layout(column-wise): [3 x n_spins]
    std::vector<float> M1(len2, 0.f);   // memory layout(column-wise): [3 x n_TE x n_spins x n_sample_length_scales]
    std::vector<float> Trajectory(len3, 0.f);   // memory layout(column-wise): [3 x n_TE x n_spins x n_sample_length_scales]

    std::cout << std::string(50, '=') << std::endl;
    for (int16_t fieldmap_no=0; fieldmap_no<param.n_fieldmaps; fieldmap_no++)
    {
        bool hasXYZ0 = false;
        // ========== load files (field-maps, xyz0, m0) ==========
        if(file_utils::read_fieldmap(filenames.at("fieldmap")[fieldmap_no], fieldmap, mask, param) == false)
            return false;

        if(filenames.at("xyz0")[fieldmap_no].empty() == false)
        {
            if(file_utils::read_file(filenames.at("xyz0")[fieldmap_no], XYZ0) == false)
                return false;
            
            if(param.enMultiTissue == false)
            {
                BOOST_LOG_TRIVIAL(info)  << "Checking XYZ0 is not in the mask..." << std::endl;
                if(is_masked(XYZ0, mask, &param))
                {
                    BOOST_LOG_TRIVIAL(error)  << "Elements of XYZ0 are in the mask or out of range. Aborting...!" << std::endl;
                    return false;
                }
            }
            hasXYZ0 = true;
        }

        if(filenames.at("m0")[fieldmap_no].empty() == false)
        {
            if(file_utils::read_file(filenames.at("m0")[fieldmap_no], M0) == false)
            {
                return false;
            }
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
            param.scale2grid[i] = (param.fieldmap_size[i] - 1.) / param.sample_length[i];
        
        if (hasXYZ0 && param.n_sample_length_scales > 1)
        {
            BOOST_LOG_TRIVIAL(error) << "loading XYZ0 from file while having more than 1 sample length scales is not supported!" << std::endl;
            return false;
        }

        // ========== distributing between devices ==========
        std::vector<float *> d_pFieldMap(device_count, NULL);
        std::vector<float *> d_M0(device_count, NULL), d_M1(device_count, NULL);
        std::vector<float *> d_XYZ1(device_count, NULL), d_XYZ0(device_count, NULL), d_XYZ0_scaled(device_count, NULL);
        std::vector<uint8_t *>  d_pMask(device_count, NULL);
        std::vector<simulation_parameters *> d_param(device_count, NULL);
        std::vector<cudaStream_t> streams(device_count, NULL);

        #pragma omp parallel for
        for(uint32_t d=0; d<device_count; d++)
        {  
            simulation_parameters param_local;
            memcpy(&param_local, &param, sizeof(simulation_parameters));
            param_local.seed += d * param.n_spins; // different seed for each GPU

            checkCudaErrors(cudaSetDevice(d));            
            checkCudaErrors(cudaStreamCreate(&streams[d]));
            // allocate memory on GPU
            checkCudaErrors(cudaMalloc((void**)&d_param[d],         sizeof(simulation_parameters)));
            checkCudaErrors(cudaMalloc((void**)&d_pFieldMap[d],     sizeof(fieldmap[0]) * fieldmap.size()));   
            checkCudaErrors(cudaMalloc((void**)&d_pMask[d],         sizeof(mask[0]) * mask.size())); 
            checkCudaErrors(cudaMalloc((void**)&d_XYZ0[d],          sizeof(float) * 3 * param.n_spins));
            checkCudaErrors(cudaMalloc((void**)&d_XYZ0_scaled[d],   sizeof(float) * 3 * param.n_spins));
            checkCudaErrors(cudaMalloc((void**)&d_XYZ1[d],          sizeof(float) * 3 * param.n_spins));
            checkCudaErrors(cudaMalloc((void**)&d_M0[d],            sizeof(float) * 3 * param.n_spins));
            checkCudaErrors(cudaMalloc((void**)&d_M1[d],            sizeof(float) * 3 * param.n_TE * param.n_spins));
            // copy data to GPU
            checkCudaErrors(cudaMemcpyAsync(d_pFieldMap[d], fieldmap.data(),        fieldmap.size()*sizeof(fieldmap[0]), cudaMemcpyHostToDevice, streams[d]));
            checkCudaErrors(cudaMemcpyAsync(d_pMask[d],     mask.data(),            mask.size() * sizeof(mask[0]),       cudaMemcpyHostToDevice, streams[d]));
            checkCudaErrors(cudaMemcpyAsync(d_param[d],     &param_local,           sizeof(simulation_parameters),       cudaMemcpyHostToDevice, streams[d]));
            checkCudaErrors(cudaMemcpyAsync(d_M0[d],        &M0[3*param.n_spins*d], 3*param.n_spins*sizeof(M0[0]),       cudaMemcpyHostToDevice, streams[d]));

            // convert fieldmap from Tesla to degree per dwell time
            cu_scaleArray<<<uint64_t(fieldmap.size()/THREADS_PER_BLOCK)+1, THREADS_PER_BLOCK, 0, streams[d]>>>(d_pFieldMap[d], param.B0*GAMMA*param.dt*RAD2DEG, fieldmap.size());
            if(hasXYZ0 == false)
            {   // generate initial spatial position for spins, based on sample_length_ref
                BOOST_LOG_TRIVIAL(info) << "GPU" << d << " Generating random initial position for spins... ";
                cu_randPosGen<<<numBlocks, THREADS_PER_BLOCK, 0, streams[d]>>>(d_XYZ0[d], d_param[d], d_pMask[d]);
                gpuCheckKernelExecutionError( __FILE__, __LINE__);
                BOOST_LOG_TRIVIAL(info) << "GPU" << d << " Done!";
            }
            else // copy initial spatial position and magnetization for spins
                checkCudaErrors(cudaMemcpyAsync(d_XYZ0[d], &XYZ0[3*param.n_spins*d], 3*param.n_spins*sizeof(XYZ0[0]), cudaMemcpyHostToDevice, streams[d]));  
            checkCudaErrors(cudaStreamSynchronize(streams[d]));    
        }

        // ========== run ==========        
        cudaEvent_t start;
        cudaEvent_t end;
        checkCudaErrors(cudaEventCreate(&start));
        checkCudaErrors(cudaEventCreate(&end));
        checkCudaErrors(cudaEventRecord(start));
        
        tqdm bar;
        simulation_parameters param_local;
        memcpy(&param_local, &param, sizeof(simulation_parameters));
        for (int32_t sl = 0; sl < param.n_sample_length_scales; sl++)
        {
            for (int i = 0; i < 3; i++)
            {
                param_local.sample_length[i] = sample_length_scales[sl] * param.sample_length[i];
                param_local.scale2grid[i] = (param_local.fieldmap_size[i] - 1.) / param_local.sample_length[i];
            }
            
            #pragma omp parallel for
            for (int32_t d = 0; d < device_count; d++)
            {                    
                checkCudaErrors(cudaSetDevice(d));
                cudaMemcpy(d_param[d], &param_local, sizeof(simulation_parameters), cudaMemcpyHostToDevice);

                cu_scalePos<<<numBlocks, THREADS_PER_BLOCK, 0, streams[d]>>>(d_XYZ0_scaled[d], d_XYZ0[d], sample_length_scales[sl], param.n_spins);
                gpuCheckKernelExecutionError(__FILE__, __LINE__);
                
                cu_sim<<<numBlocks, THREADS_PER_BLOCK, 0, streams[d]>>>(d_param[d], d_pFieldMap[d], d_pMask[d], d_M0[d], d_XYZ0_scaled[d], d_M1[d], d_XYZ1[d]);
                gpuCheckKernelExecutionError(__FILE__, __LINE__);

                int shift = 3*param.n_TE*param.n_spins*device_count*sl + 3*param.n_TE*param.n_spins*d;
                checkCudaErrors(cudaMemcpyAsync(M1.data()   + shift, d_M1[d]  , 3*param.n_TE*param.n_spins*sizeof(float), cudaMemcpyDeviceToHost, streams[d]));
                shift = 3*param.n_spins*device_count*sl + 3*param.n_spins*d;
                checkCudaErrors(cudaMemcpyAsync(XYZ1.data() + shift, d_XYZ1[d], 3*param.n_spins*sizeof(float), cudaMemcpyDeviceToHost, streams[d]));
            }
            bar.progress(sl, param.n_sample_length_scales);
        }
        bar.finish();

        float elapsedTime;
        checkCudaErrors(cudaEventRecord(end));
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, end));
        std::cout << "Entire simulation over " << device_count << " GPU(s) took " << std::fixed << std::setprecision(2) << elapsedTime/1000. << " second(s)" << std::endl;

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
        
        // ========== save results ========== 
        file_utils::output_header hdr(3, param.n_TE, param.n_spins * device_count, param.n_sample_length_scales);
        file_utils::save_output(M1, filenames.at("m1")[fieldmap_no], hdr, sample_length_scales);

        hdr.dim2 = 1;
        if(filenames.at("xyz1")[fieldmap_no].empty() == false) // do not save if filename is empty
            file_utils::save_output(XYZ1, filenames.at("xyz1")[fieldmap_no], hdr, sample_length_scales);

        std::cout << std::string(50, '=') << std::endl;
    }
    return true;
} 


bool dump_settings(simulation_parameters param, map<string, vector<string> > filenames, std::vector<float> sample_length_scales)
{
    std::stringstream ss;
    ss  << "Dumping settings:" << std::endl;
    for (std::map<std::string, std::vector<std::string>>::iterator it=filenames.begin(); it!=filenames.end(); ++it)
        for (int i = 0; i< it->second.size(); i++)
            ss << it->first << "[" << i << "] = " << it->second.at(i) << std::endl;
    
    ss<< "\nSample length scale = [";
    for (int32_t i = 0; i < param.n_sample_length_scales; i++)
        ss << sample_length_scales[i] << ", ";
    ss << "\b\b]\n" << std::endl;

    file_utils::input_header hdr_in;
    if(file_utils::read_header(filenames.at("fieldmap")[0], hdr_in) == false)
    {
        BOOST_LOG_TRIVIAL(error) << "reading header of fieldmap " << filenames.at("fieldmap")[0] << " failed. Aborting...!";
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
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
    boost::program_options::notify(vm);

    boost::log::add_file_log(boost::log::keywords::file_name=LOG_FILE, boost::log::keywords::target_file_name = LOG_FILE, boost::log::keywords::format = "[%TimeStamp%] [%Severity%]: %Message%");
    boost::log::add_common_attributes();
    
    // ========== print help ==========
    if (vm.count("help") || vm.count("configs") == 0 || argc == 1)
    {
        std::cout << desc;
        print_device_info();
        return 1;
    }

    std::vector<std::string> config_files = vm["configs"].as<std::vector<std::string>>();  
    bool bNoSim = vm.count("sim_off") > 0;

    // ========== loop over configs and simulate ==========
    std::cout << "Running simulation for " << config_files.size() << " config(s)..." << std::endl;
    for(const auto& cfile : config_files)
    {
        map<string, vector<string> > filenames = {{"fieldmap", 	vector<string>()},  // input:  map of off-resonance in Tesla
                                                  {"xyz0", 		vector<string>()},  // input:  spins starting spatial positions in meters
                                                  {"xyz1", 		vector<string>()},  // output: spins last spatial positions in meters
                                                  {"m0", 		vector<string>()},  // input:  spins initial magnetization
                                                  {"m1", 		vector<string>()}}; // output: spins final magnetization

        std::vector<float> sample_length_scales;
        simulation_parameters param;

        // ========== read config file ==========
        param.fieldmap_size[0] = param.fieldmap_size[1] = param.fieldmap_size[2] = 0;
        param.sample_length[0] = param.sample_length[1] = param.sample_length[2] = 0.f;
        bStatus &= file_utils::read_config(cfile, param, sample_length_scales, filenames);
 
        if (param.seed == 0)
            param.seed = std::random_device{}();

        param.n_timepoints = param.TR / param.dt; // includes start point

        // ========== dump settings ==========
        bStatus &= dump_settings(param, filenames, sample_length_scales);
        // ========== GPU memory is enough ==========
        bStatus &= check_memory_size(param.get_required_memory(getDeviceCount()));
        if (bStatus == false)
        {
            std::cout << ERR_MSG << "Simulation failed. See the log file: " << LOG_FILE <<". Aborting...!" << std::endl;
            return 1;
        }
        
        if (bNoSim)
            continue;

        if(simulate(param, filenames, sample_length_scales) == false)
        {
            std::cout << ERR_MSG << "Simulation failed. See the log file: " << LOG_FILE <<". Aborting...!" << std::endl;
            return 1;
        }
    }
    std::cout << "Simulation(s) finished successfully!" << std::endl;
    return 0;
}
