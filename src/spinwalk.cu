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
#include <algorithm>
#include <execution>
#include <filesystem>
#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include "indicators.hpp"
#include "version.h"
#include "kernels.cuh"
#include "file_utils.h"
#include "shapes/cylinder.cuh"
#include "shapes/sphere.cuh"

#ifdef __CUDACC__
#include <cuda_runtime.h>
#include "helper_cuda.h"
#include "helper.cuh"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/iterator/constant_iterator.h>
#endif

#define THREADS_PER_BLOCK  64

namespace bl = boost::log;
using namespace indicators;

bool run(simulation_parameters param, std::map<std::string, std::vector<std::string> > filenames, std::vector<double> fov_scale)
{
    auto start_config = std::chrono::high_resolution_clock::now();
    int64_t old_elapsed = 0;
    // ========== checking number of GPU(s) ==========
    int32_t device_count=1;
#ifdef __CUDACC__
    checkCudaErrors(cudaGetDeviceCount(&device_count));
    BOOST_LOG_TRIVIAL(info) << "Number of available GPU(s): " << device_count; 
#endif    
    // param.n_spins /= device_count; // spins will be distributed in multiple GPUs (if there is). We suppose it is divisible 
    size_t numBlocks = (param.n_spins + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    
    // ========== allocate memory on CPU ==========
    size_t trj  = param.enRecordTrajectory ? param.n_timepoints * (param.n_dummy_scan + 1) : 1;
    size_t len0 = 3 * param.n_spins;
    size_t len1 = len0 * param.n_fov_scale * trj;
    size_t len2 = len0 * param.n_fov_scale * param.n_TE;
    BOOST_LOG_TRIVIAL(info) << "Memory size: XYZ0=" << len0 << ", XYZ1=" << len1 << ", M0=" << len0 << ", M1=" << len2 << "";
    std::vector<float>      fieldmap(param.fieldmap_exist?param.matrix_length:0, 0.f);
    std::vector<uint8_t>    mask(param.matrix_length, 0);
    std::vector<float>      XYZ0(len0, 0.f);     // memory layout(row-major): [n_spins x 3]
    std::vector<float>      XYZ0_scaled(len0, 0.f);     // memory layout(row-major): [n_spins x 3]
    std::vector<float>      XYZ1(len1, 0.f);     // memory layout(row-major): [n_fov_scale x n_spins x timepoints x 3] or [n_fov_scale x n_spins x 1 x 3]
    std::vector<float>      M0(len0, 0.f);       // memory layout(row-major): [n_spins x 3]
    std::vector<float>      M1(len2, 0.f);       // memory layout(row-major): [n_fov_scale x n_spins x n_TE x 3]
    std::vector<uint8_t>    T(M1.size()/3, 0);   // memory layout(row-major): [n_fov_scale x n_spins x n_TE x 1]
    
    // ========== allocate memory on GPU ==========
#ifdef __CUDACC__
    thrust::device_vector<float>   d_pFieldMap;
    thrust::device_vector<float>   d_M0(M0.size());
    thrust::device_vector<float>   d_M1(M1.size() / param.n_fov_scale);
    thrust::device_vector<float>   d_XYZ1(XYZ1.size() / param.n_fov_scale); 
    thrust::device_vector<float>   d_XYZ0(XYZ0.size());
    thrust::device_vector<float>   d_XYZ0_scaled(XYZ0_scaled.size());
    thrust::device_vector<uint8_t> d_T(T.size() / param.n_fov_scale);
    thrust::device_vector<uint8_t> d_pMask(mask.size());
    simulation_parameters *d_param = nullptr;
    checkCudaErrors(cudaMalloc(&d_param, sizeof(simulation_parameters)));
    if(param.fieldmap_exist)
        d_pFieldMap.resize(fieldmap.size());
#endif

    for (int16_t fieldmap_no=0; fieldmap_no<param.n_fieldmaps; fieldmap_no++)
    {
        // ========== load files (field-maps, xyz0, m0) ==========
        std::cout << "Loading phantom: " << std::filesystem::path(filenames.at("phantom")[fieldmap_no]).filename().string() << "\n";
        if(file_utils::read_phantom(filenames.at("phantom")[fieldmap_no], fieldmap, mask, &param) == false)
            return false;
        param.matrix_length = mask.size(); // update the matrix length based on the mask size from the recent read
#ifdef __CUDACC__
        if(d_pMask.size() != mask.size() && param.no_gpu == false)
            d_pMask.resize(mask.size());
        if(d_pFieldMap.size() != fieldmap.size() && param.fieldmap_exist && param.no_gpu == false)
            d_pFieldMap.resize(fieldmap.size());
#endif
        // convert fieldmap from T to degree per timestep
        float Tesla2deg = param.B0 * param.timestep_us * 1e-6 * GAMMA * RAD2DEG;
        if(param.fieldmap_exist) 
            std::transform(std::execution::par_unseq, fieldmap.begin(), fieldmap.end(), fieldmap.begin(), [Tesla2deg](float x) { return x * Tesla2deg;});

        if(fieldmap_no < filenames.at("xyz0").size())
        {   
            uint32_t n_spins = product(file_utils::get_size_h5(filenames.at("xyz0")[fieldmap_no], "/XYZ")) / 3;            
            if(n_spins != param.n_spins)
            {
                BOOST_LOG_TRIVIAL(error) << "Number of spins in XYZ0 file does not match with the number of spins in the config file! " << n_spins << " vs " << param.n_spins ;
                return false;
            }
            if(file_utils::read_h5(filenames.at("xyz0")[fieldmap_no], XYZ0.data(), "/XYZ", "float") == false)
                return false;
        }
        else
        {   // generate initial spatial position for spins
            BOOST_LOG_TRIVIAL(info) << "GPU " << 1 << ") Generating random initial position for spins... (seed = " << param.seed << ", GPU grid = " << numBlocks << " x " << THREADS_PER_BLOCK << ")";
            randPosGen(XYZ0.data(), param);
            BOOST_LOG_TRIVIAL(info) << "GPU " << 1 << ") Done!";
        }

        if(fieldmap_no < filenames.at("m0").size())
        {
            uint32_t n_spins = product(file_utils::get_size_h5(filenames.at("m0")[fieldmap_no], "/M")) / 3; 
            if(n_spins != param.n_spins)
            {
                BOOST_LOG_TRIVIAL(error) << "Number of spins in M0 file does not match with the number of spins in the config file! " << n_spins << " vs " << param.n_spins ;
                return false;
            }
            if(file_utils::read_h5(filenames.at("m0")[fieldmap_no], M0.data(), "/M", "float") == false)
                return false;
        }
        else
        {   // all spins are aligned with B0 (M0 = (0, 0, 1))
            uint32_t index = 0;
            BOOST_LOG_TRIVIAL(info) << "Generating M0(0, 0, 1)..." << std::endl;
            std::generate(M0.begin(), M0.end(), [&index](){return (index++ % 3 == 2) ? 1.f : 0.f;});
        }

        for(int i=0; i<M0.size()/3; i += M0.size()/3/2)
            BOOST_LOG_TRIVIAL(info) << "M0 of the spin " << i << " = (" << M0[3*i] << ", " << M0[3*i+1] << ", " << M0[3*i+2] << ")" << std::endl;
        
        for(int i=0; i<3; i++)
            param.scale2grid[i] = param.fieldmap_size[i] / param.fov[i];
        
#ifdef __CUDACC__
        // ========== copy data to GPU(s) ==========      
        if (param.no_gpu == false)
        {
            checkCudaErrors(cudaMemcpyAsync(d_param, &param, sizeof(simulation_parameters), cudaMemcpyHostToDevice));
            d_pMask = mask; 
            d_M0    = M0;
            d_XYZ0  = XYZ0;           
            if(param.fieldmap_exist)    
                d_pFieldMap = fieldmap;
        }
#endif
        // ========== run ==========   
        auto start_sim = std::chrono::high_resolution_clock::now();
        ProgressBar bar{option::ShowPercentage{true}, option::Start{"["}, option::Fill{"="}, option::Lead{">"}, option::End{"]"}};
        simulation_parameters param_local;
        memcpy(&param_local, &param, sizeof(simulation_parameters));
        std::vector<uint32_t> v(param_local.n_spins);
        std::generate(std::execution::seq, v.begin(), v.end(), [n = 0] () mutable { return n++; });  
        for (int32_t sl = 0; sl < param.n_fov_scale; sl++)
        {
            // virtual FoV scaling
            for (int i = 0; i < 3; i++)
            {
                param_local.fov[i] = fov_scale[sl] * param.fov[i];
                param_local.scale2grid[i]    = param_local.fieldmap_size[i] / param_local.fov[i];
                for (int n = 0; n < param.n_tissue_type; n++)
                {
                    double min_convergence = 0.95 * param.diffusivity[n] * sqrt(1.0 * param.TE_us[param.n_TE-2]); // https://submissions.mirasmart.com/ISMRM2024/Itinerary/PresentationDetail.aspx?evdid=4684
                    if(param_local.fov[i] < min_convergence )
                    {
                        BOOST_LOG_TRIVIAL(error) << "Virtual FoV (= " << param_local.fov[i] << ") is smaller than minimum convergence length (= " << min_convergence << ")!";
                        BOOST_LOG_TRIVIAL(error) << "Original FoV (= " << param.fov[i] << ") FoV Scale (= " << fov_scale[sl]  << ")!";
                        return false;
                    }
                }
            }           
            // run simulation kernel
            if(param.no_gpu)
            {
                float scale = fov_scale[sl];                
                // scale position to mimic the different volume size
                std::transform(std::execution::par_unseq, XYZ0.begin(), XYZ0.end(), XYZ0_scaled.begin(), [scale](auto& c){return c*scale;}); 
                std::for_each(std::execution::par_unseq, v.begin(), v.end(), [&](int spin) {sim(&param_local, fieldmap.data(), 
                                                                                                              mask.data(), 
                                                                                                              M0.data(),
                                                                                                              XYZ0_scaled.data(), 
                                                                                                              M1.data()   + 3*param.n_TE*param.n_spins*sl, 
                                                                                                              XYZ1.data() + 3*param.n_spins*trj*sl, 
                                                                                                              T.data()    + param.n_TE*param.n_spins*sl,
                                                                                                              spin);});
            }
#ifdef __CUDACC__  
            else
            {
            BOOST_LOG_TRIVIAL(info) << "GPU " << 1 << ") Fieldmap " << fieldmap_no << ", simulating sample length scale " << fov_scale[sl];
            checkCudaErrors(cudaMemcpy(d_param, &param_local, sizeof(simulation_parameters), cudaMemcpyHostToDevice));
            // scale position to mimic the different volume size
            thrust::transform(d_XYZ0.begin(), d_XYZ0.end(), thrust::make_constant_iterator(fov_scale[sl]), d_XYZ0_scaled.begin(), thrust::multiplies<float>());
            gpuCheckKernelExecutionError(__FILE__, __LINE__);            
            cu_sim<<<numBlocks, THREADS_PER_BLOCK, 0>>>(d_param,thrust::raw_pointer_cast(d_pFieldMap.data()), 
                                                                thrust::raw_pointer_cast(d_pMask.data()), 
                                                                thrust::raw_pointer_cast(d_M0.data()), 
                                                                thrust::raw_pointer_cast(d_XYZ0_scaled.data()), 
                                                                thrust::raw_pointer_cast(d_M1.data()), 
                                                                thrust::raw_pointer_cast(d_XYZ1.data()), 
                                                                thrust::raw_pointer_cast(d_T.data())); 
            gpuCheckKernelExecutionError(__FILE__, __LINE__);
            // copy data back to CPU
            size_t shift = 3*param.n_TE*param.n_spins*sl;
            thrust::copy(d_M1.begin(), d_M1.end(), M1.begin() + shift);
            shift = 3*param.n_spins*trj*sl;
            thrust::copy(d_XYZ1.begin(), d_XYZ1.end(), XYZ1.begin() + shift);
            shift = param.n_TE*param.n_spins*sl;
            thrust::copy(d_T.begin(), d_T.end(), T.begin() + shift);
            }
#endif
            // bar.progress(sl, param.n_fov_scale);
            bar.set_progress(100 * (sl+1)/float(param.n_fov_scale));
        }
        
        auto end_config     = std::chrono::high_resolution_clock::now();        
        auto elapsed_sim    = std::chrono::duration_cast<std::chrono::milliseconds>(end_config - start_sim).count() / 1000.0;
        auto elapsed_config = std::chrono::duration_cast<std::chrono::seconds>(end_config - start_config).count();
        int precision = elapsed_sim>10 ? 0 : (elapsed_sim > 1 ? 1 : 3);
        std::cout << "Simulation took " << std::fixed << std::setprecision(precision) <<  elapsed_sim << " sec., everything else took " << elapsed_config - elapsed_sim - old_elapsed<< " sec.\n";
        old_elapsed = elapsed_config;

        // ========== save results ========== 
        std::cout << "Saving the results to disk." << "\n";
        std::string f = filenames.at("output")[fieldmap_no];
        if (std::filesystem::exists(f)) 
        {
            std::filesystem::remove(f);
            BOOST_LOG_TRIVIAL(info) << "File " << f << " already exists. Overwriting it.";
        }
        
        std::vector<size_t> dims = {param.n_fov_scale, param.n_spins, param.n_TE, 3};
        file_utils::save_h5(f, M1.data(), dims, "M", "float");

        dims[2] = param.enRecordTrajectory ? param.n_timepoints * (param.n_dummy_scan + 1) : 1;
        file_utils::save_h5(f, XYZ1.data(), dims, "XYZ", "float");

        dims[3] = 1; dims[2] = param.n_TE;
        file_utils::save_h5(f, T.data(), dims, "T", "uint8_t");

        dims[0] = fov_scale.size(); dims[1] = 1; dims[2] = 1; dims[3] = 1;
        file_utils::save_h5(f, fov_scale.data(), dims, "scales", "double");
    }

#ifdef __CUDACC__ 
    // ============ clean up GPU ============
    checkCudaErrors(cudaFree(d_param));        
#endif

    return true;
} 


bool dump_settings(simulation_parameters param, std::map<std::string, std::vector<std::string> > filenames, std::vector<double> fov_scale)
{
    std::stringstream ss;
    ss << "Dumping settings:" << '\n';
    for (std::map<std::string, std::vector<std::string>>::iterator it=filenames.begin(); it!=filenames.end(); ++it)
        for (int i = 0; i< it->second.size(); i++)
            ss << it->first << "[" << i << "] = " << it->second.at(i) << '\n';
    
    ss << "\nFoV scale = [";
    for (int32_t i = 0; i < param.n_fov_scale; i++)
        ss << fov_scale[i] << ", ";
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
    int32_t phantom_resolution, phantom_seed, device_id = 0;
    std::vector<std::string>  config_files;
    print_logo();

    // ========== parse command line arguments ==========
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
        ("help,h", "help message (this menu)")
        ("configs,c", po::value<std::vector<std::string>>(&config_files)->multitoken(), "config. files as many as you want. e.g. -c config1.ini config2.ini ... configN.ini")
#ifdef __CUDACC__
        ("use_cpu,p", "only run on CPU (default: GPU)")
        ("device,g",po::value<int32_t>(&device_id)->default_value(0), "select GPU device (if there are multiple GPUs)")
#endif
        ("cylinder,C", "generate phantom filled with cylinders")
        ("sphere,S", "generate phantom filled with spheres")
        ("orientation,o", po::value<float>(&phantom_orientation)->default_value(90.0), "orientation of the cylinders in degree with respect to B0")
        ("radius,r", po::value<float>(&phantom_radius)->default_value(50), "radius of the cylinders/spheres in um (negative value = random radius)")
        ("seed,s", po::value<int32_t>(&phantom_seed)->default_value(-1), "seed for random number generator in phantom creator (negative value = random seed)")
        ("BVF,b", po::value<float>(&phantom_BVF)->default_value(10.0), "fraction of shapes to entire volume <0.0 100.0> (i.e. blood volume fraction)")
        ("fov,v", po::value<float>(&phantom_fov)->default_value(1000.0), "voxel field of view in um (isotropic)")
        ("resolution,z", po::value<int32_t>(&phantom_resolution)->default_value(500), "base resolution")
        ("dchi,d", po::value<float>(&phantom_dchi)->default_value(0.11e-6), "susceptibility difference between fully deoxygenated blood (inside cylinders/spheres) and tissue (outside cylinders/spheres) (default: 0.11e-6 in cgs units)")
        ("oxy_level,y", po::value<float>(&phantom_oxy_level)->default_value(0.75), "blood oxygenetation level <0.0 1.0> (negative value = exclude off-resonance effect and only generate the mask)")
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
    
    // ========== setup log ==========
    std::string log_filename = "spinwalk_" + std::to_string(device_id) + ".log";
    auto fileSink = bl::add_file_log(bl::keywords::file_name=log_filename, bl::keywords::target_file_name = log_filename, bl::keywords::format = "[%TimeStamp%] [%Severity%]: %Message%", bl::keywords::auto_flush = true);
    bl::add_common_attributes();
    
    // ========== print help ==========
    if (vm.count("help") || argc == 1 || unreg.size() > 0)
    {
        std::cout << desc;
        print_device_info();
        return 0;
    }

    std::cout << "Log file location: " << std::filesystem::current_path() / log_filename << '\n';

    // ========== generate phantom ==========
    if (vm.count("cylinder"))
    {
        cylinder cyl(phantom_fov, phantom_resolution, phantom_dchi, phantom_oxy_level, phantom_radius, phantom_BVF, phantom_orientation, phantom_seed, phantom_output);
        cyl.run();
    }
    if (vm.count("sphere"))
    {
        sphere sph(phantom_fov, phantom_resolution, phantom_dchi, phantom_oxy_level, phantom_radius, phantom_BVF, phantom_seed, phantom_output);
        sph.run();
    }

    // ========== loop over configs and simulate ==========
    if (config_files.size() == 0)
        return 0;
    
    std::cout << "Running simulation for " << config_files.size() << " config(s)..." << "\n\n";
    auto start = std::chrono::high_resolution_clock::now();
    for(const auto& cfile : config_files)
    {
        std::cout << "<" << std::filesystem::path(cfile).filename().string() << ">\n";
        std::map<std::string, std::vector<std::string> > filenames = {{"phantom", 	std::vector<std::string>()},  // input:  map of off-resonance in Tesla
                                                                      {"xyz0", 		std::vector<std::string>()},  // input:  spins starting spatial positions in meters
                                                                      {"m0", 		std::vector<std::string>()},  // input:  spins initial magnetization
                                                                      {"output", 	std::vector<std::string>()},  // output: spins final magnetization + spatial positions in meters + tissue index
                                                                     };

        std::vector<double> fov_scale;
        simulation_parameters param;
        param.no_gpu = vm.count("use_cpu") > 0;
        
        // ========== read config file ==========
        bStatus &= file_utils::read_config(cfile, &param, fov_scale, filenames);
            
        // ========== dump settings ==========
        bStatus &= dump_settings(param, filenames, fov_scale);

        // ========== Check GPU memory ==========
        bStatus &= check_memory_size(param.get_required_memory(getDeviceCount()));

        if (bStatus == false)
        {
            std::cout << ERR_MSG << "Simulation failed. See the log file " << log_filename <<", Aborting...!" << std::endl;
            return 1;
        }
        
#ifdef __CUDACC__
        if(param.no_gpu == false)
        {
            if (device_id >= getDeviceCount())
            {
                std::cout << ERR_MSG << "Device ID " << device_id << " is not available! Aborting...!" << std::endl;
                return 1;
            }
            cudaSetDevice(device_id);
        }
#endif

        std::cout << "Simulation starts..." << std::endl;
        if(run(param, filenames, fov_scale) == false)
        {
            std::cout << ERR_MSG << "Simulation failed. See the log file " << log_filename <<", Aborting...!" << std::endl;
            return 1;
        }        
        std::cout << "\n";
    }
    std::cout << "Simulation(s) finished successfully! Total elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000. << " second(s)."  << std::endl;
    return 0;
}
