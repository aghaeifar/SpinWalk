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
#include <iomanip>
#include <filesystem>
#include <cuda_runtime.h>
#include "tqdm.h"
#include "helper.cuh"
#include "kernels.cuh"
#include "file_utils.h"
#include "shapes/cylinder.cuh"
#include "shapes/sphere.cuh"
#include "helper_cuda.h"
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
    std::vector<float>      fieldmap(param.fieldmap_exist?param.matrix_length:0, 0.f);
    std::vector<uint8_t>    mask(param.matrix_length, 0);
    std::vector<float>      XYZ0(len0, 0.f);     // memory layout(column-wise): [3 x n_spins]
    std::vector<float>      XYZ1(len1, 0.f);     // memory layout(column-wise): [3 x timepoints x n_spins x n_sample_length_scales] or [3 x 1 x n_spins x n_sample_length_scales]
    std::vector<float>      M0(len0, 0.f);       // memory layout(column-wise): [3 x n_spins]
    std::vector<float>      M1(len2, 0.f);       // memory layout(column-wise): [3 x n_TE x n_spins x n_sample_length_scales]
    std::vector<uint8_t>    T(M1.size()/3, 0);   // memory layout(column-wise): [1 x n_TE x n_spins x n_sample_length_scales]
    BOOST_LOG_TRIVIAL(info) << "Memory allocation (CPU) took " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - st).count() << " ms";
    
    // ========== allocate memory on GPU ==========
    std::vector<float *>    d_pFieldMap(device_count, nullptr);
    std::vector<float *>    d_M0(device_count, nullptr), d_M1(device_count, nullptr);
    std::vector<float *>    d_XYZ1(device_count, nullptr), d_XYZ0(device_count, nullptr), d_XYZ0_scaled(device_count, nullptr);
    std::vector<uint8_t *>  d_T(device_count, nullptr);
    std::vector<uint8_t *>  d_pMask(device_count, nullptr);
    std::vector<simulation_parameters *> d_param(device_count, nullptr);
    std::vector<cudaStream_t> streams(device_count, nullptr);
    
    #pragma omp parallel for
    for(int32_t d=0; d<device_count; d++)
    {  
        checkCudaErrors(cudaSetDevice(d));            
        checkCudaErrors(cudaStreamCreate(&streams[d]));
        // allocate memory on GPU
        checkCudaErrors(cudaMalloc((void**)&d_param[d],         sizeof(simulation_parameters)));         
        checkCudaErrors(cudaMalloc((void**)&d_XYZ0[d],          sizeof(XYZ0[0]) * XYZ0.size() / device_count));
        checkCudaErrors(cudaMalloc((void**)&d_XYZ0_scaled[d],   sizeof(XYZ0[0]) * XYZ0.size() / device_count));
        checkCudaErrors(cudaMalloc((void**)&d_XYZ1[d],          sizeof(XYZ1[0]) * XYZ1.size() / device_count / param.n_sample_length_scales));
        checkCudaErrors(cudaMalloc((void**)&d_M0[d],            sizeof(M0[0]) * M0.size() / device_count));
        checkCudaErrors(cudaMalloc((void**)&d_M1[d],            sizeof(M1[0]) * M1.size() / device_count / param.n_sample_length_scales));
        checkCudaErrors(cudaMalloc((void**)&d_T[d],             sizeof(T[0])  * T.size()  / device_count / param.n_sample_length_scales));
        checkCudaErrors(cudaMalloc((void**)&d_pMask[d],         sizeof(mask[0]) * mask.size())); 
        if(param.fieldmap_exist)
            checkCudaErrors(cudaMalloc((void**)&d_pFieldMap[d], sizeof(fieldmap[0]) * fieldmap.size())); 
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
        if(file_utils::read_fieldmap(filenames.at("FIELDMAP")[fieldmap_no], fieldmap, mask, &param) == false)
            return false;

        if(fieldmap_no < filenames.at("XYZ0").size())
        {   
            XYZ0.resize(product(file_utils::get_size_h5(filenames.at("XYZ0")[fieldmap_no], "/XYZ")));            
            if(XYZ0.size()/3 != param.n_spins*device_count)
            {
                BOOST_LOG_TRIVIAL(error) << "Number of spins in XYZ0 file does not match with the number of spins in the config file! " << XYZ0.size()/3 << " vs " << param.n_spins*device_count ;
                return false;
            }
            if(file_utils::read_h5(filenames.at("XYZ0")[fieldmap_no], XYZ0.data(), "/XYZ", "float") == false)
                return false;
            hasXYZ0 = true;
        }

        if(fieldmap_no < filenames.at("M0").size())
        {
            M0.resize(product(file_utils::get_size_h5(filenames.at("M0")[fieldmap_no], "/M"))); 
            if(M0.size() != param.n_spins*device_count)
            {
                BOOST_LOG_TRIVIAL(error) << "Number of spins in M0 file does not match with the number of spins in the config file! " << M0.size() << " vs " << param.n_spins*device_count ;
                return false;
            }
            if(file_utils::read_h5(filenames.at("M0")[fieldmap_no], M0.data(), "/M", "float") == false)
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
            checkCudaErrors(cudaMemcpyAsync(d_pMask[d],     mask.data(),            mask.size() * sizeof(mask[0]),       cudaMemcpyHostToDevice, streams[d]));
            checkCudaErrors(cudaMemcpyAsync(d_param[d],     &param_local,           sizeof(simulation_parameters),       cudaMemcpyHostToDevice, streams[d]));
            checkCudaErrors(cudaMemcpyAsync(d_M0[d],        &M0[3*param.n_spins*d], 3*param.n_spins*sizeof(M0[0]),       cudaMemcpyHostToDevice, streams[d]));
            if(param.fieldmap_exist)   
            {      
                checkCudaErrors(cudaMemcpyAsync(d_pFieldMap[d], fieldmap.data(),    fieldmap.size()*sizeof(fieldmap[0]), cudaMemcpyHostToDevice, streams[d]));
                // convert fieldmap from uT to degree per timestep
                cu_scaleArray<<<uint64_t(fieldmap.size()/THREADS_PER_BLOCK)+1, THREADS_PER_BLOCK, 0, streams[d]>>>(d_pFieldMap[d], param.B0 * param.timestep_us * 1e-6 * GAMMA * RAD2DEG, fieldmap.size());
            }
            if(hasXYZ0 == false)
            {   // generate initial spatial position for spins, based on sample_length_ref
                BOOST_LOG_TRIVIAL(info) << "GPU " << d << ") Generating random initial position for spins... (seed = " << param_local.seed << ", GPU grid = " << numBlocks << " x " << THREADS_PER_BLOCK << ")";
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
        std::string f = filenames.at("output")[fieldmap_no];
        
        std::vector<size_t> dims = {3, param.n_TE, param.n_spins * device_count, param.n_sample_length_scales};
        file_utils::save_h5(f, M1.data(), dims, "M", "float");

        dims[1] = param.enRecordTrajectory ? param.n_timepoints * (param.n_dummy_scan + 1) : 1;
        file_utils::save_h5(f, XYZ1.data(), dims, "XYZ", "float");

        dims[0] = 1; dims[1] = param.n_TE;
        file_utils::save_h5(f, T.data(), dims, "T", "uint8_t");

        dims[0] = sample_length_scales.size(); dims[1] = 1; dims[2] = 1; dims[3] = 1;
        file_utils::save_h5(f, sample_length_scales.data(), dims, "scales", "double");

        std::cout << std::string(50, '=') << std::endl;
    }

    // ========== clean up GPU ==========
    #pragma omp parallel for
    for(int32_t d=0; d<device_count; d++)
    {
        checkCudaErrors(cudaSetDevice(d));   
        checkCudaErrors(cudaFree(d_param[d]));        
        checkCudaErrors(cudaFree(d_XYZ0[d]));
        checkCudaErrors(cudaFree(d_XYZ0_scaled[d]));
        checkCudaErrors(cudaFree(d_M1[d]));
        checkCudaErrors(cudaFree(d_XYZ1[d]));
        checkCudaErrors(cudaFree(d_T[d]));
        checkCudaErrors(cudaFree(d_pMask[d]));
        if(param.fieldmap_exist) 
            checkCudaErrors(cudaFree(d_pFieldMap[d]));
        checkCudaErrors(cudaStreamDestroy(streams[d]));            
    }
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(end));
    
    return true;
} 


bool dump_settings(simulation_parameters param, std::map<std::string, std::vector<std::string> > filenames, std::vector<double> sample_length_scales)
{
    std::stringstream ss;
    ss << "Dumping settings:" << '\n';
    for (std::map<std::string, std::vector<std::string>>::iterator it=filenames.begin(); it!=filenames.end(); ++it)
        for (int i = 0; i< it->second.size(); i++)
            ss << it->first << "[" << i << "] = " << it->second.at(i) << '\n';
    
    ss << "\nSample length scale = [";
    for (int32_t i = 0; i < param.n_sample_length_scales; i++)
        ss << sample_length_scales[i] << ", ";
    ss << "]\n";
    ss << param.dump();
    BOOST_LOG_TRIVIAL(info) << ss.str();
    return true;
}


int main(int argc, char * argv[])
{
    bool bStatus = true;
    std::string phantom_output;
    float phantom_radius, phantom_fov, phantom_dchi, phantom_oxy_level, phantom_orientation, phantom_BVF;
    int32_t phantom_resolution;
    std::vector<std::string>  config_files;
    print_logo();
    // ========== parse command line arguments ==========
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
        ("help,h", "help message (this menu)")
        ("configs,c", po::value<std::vector<std::string>>(&config_files)->multitoken(), "config. files as many as you want. e.g. -c config1.ini config2.ini ... configN.ini")
        ("cylinder,l", "generate phantom filled with cylinders")
        ("sphere,s", "generate phantom filled with spheres")
        ("orientation,o", po::value<float>(&phantom_orientation)->default_value(-1.0), "orientation of the cylinders in degree with respect to B0 (negative value = random orientation)")
        ("radius,r", po::value<float>(&phantom_radius)->default_value(50), "radius of the cylinders/spheres in um (negative value = random radius)")
        ("blood_volume,b", po::value<float>(&phantom_BVF)->default_value(10.0), "fraction of shapes to entire volume <0.0 100.0> ")
        ("fov,f", po::value<float>(&phantom_fov)->default_value(1000.0), "voxel field of view in um (isotropic)")
        ("resolution,z", po::value<int32_t>(&phantom_resolution)->default_value(500), "base resolution")
        ("dchi,d", po::value<float>(&phantom_dchi)->default_value(0.11e-6), "susceptibility difference between fully deoxygenated blood (inside cylinders/spheres) and tissue (outside cylinders/spheres) (default: 0.11e-6 in cgs units)")
        ("oxy_level,Y", po::value<float>(&phantom_oxy_level)->default_value(0.75), "blood oxygenetation level <0.0 1.0> (negative value = exclude off-resonance effect)")
        ("output_file,f", po::value<std::string>(&phantom_output)->default_value("./phantom.h5"), "path to save phantom (h5 format)");

    po::variables_map vm;
    std::vector<std::string> unreg;
    try{ 
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
        po::store(parsed, vm);   
        po::notify(vm); // sets default values
        unreg = po::collect_unrecognized(parsed.options, po::include_positional);
    }catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    } 
    
    auto fileSink = bl::add_file_log(bl::keywords::file_name=LOG_FILE, bl::keywords::target_file_name = LOG_FILE, bl::keywords::format = "[%TimeStamp%] [%Severity%]: %Message%", bl::keywords::auto_flush = true);
    bl::add_common_attributes();
    
    // ========== print help ==========
    if (vm.count("help") || argc == 1 || unreg.size() > 0)
    {
        std::cout << desc;
        print_device_info();
        return 0;
    }

    std::cout << "Log file location: " << std::filesystem::current_path() / LOG_FILE << '\n';

    // ========== generate phantom ==========
    if (vm.count("cylinder"))
    {
        cylinder cyl(phantom_fov, phantom_resolution, phantom_dchi, phantom_oxy_level, phantom_radius, phantom_BVF, phantom_orientation, phantom_output);
        cyl.run();
    }
    if (vm.count("sphere"))
    {
        sphere sph(phantom_fov, phantom_resolution, phantom_dchi, phantom_oxy_level, phantom_radius, phantom_BVF, phantom_output);
        sph.run();
    }

    // ========== loop over configs and simulate ==========
    if (config_files.size() == 0)
        return 0;
    
    std::cout << "Running simulation for " << config_files.size() << " config(s)..." << '\n';
    auto start = std::chrono::steady_clock::now();
    for(const auto& cfile : config_files)
    {
        std::map<std::string, std::vector<std::string> > filenames = {{"FIELDMAP", 	std::vector<std::string>()},  // input:  map of off-resonance in Tesla
                                                                      {"XYZ0", 		std::vector<std::string>()},  // input:  spins starting spatial positions in meters
                                                                      {"M0", 		std::vector<std::string>()},  // input:  spins initial magnetization
                                                                      {"output", 	std::vector<std::string>()},  // output: spins final magnetization + spatial positions in meters + tissue index
                                                                     };

        std::vector<double> sample_length_scales;
        simulation_parameters param;
        
        // ========== read config file ==========
        bStatus &= file_utils::read_config(cfile, &param, sample_length_scales, filenames);
            
        // ========== dump settings ==========
        bStatus &= dump_settings(param, filenames, sample_length_scales);

        // ========== Check GPU memory ==========
        bStatus &= check_memory_size(param.get_required_memory(getDeviceCount()));

        if (bStatus == false)
        {
            std::cout << ERR_MSG << "Simulation failed. See the log file " << LOG_FILE <<", Aborting...!" << std::endl;
            return 1;
        }
        
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
