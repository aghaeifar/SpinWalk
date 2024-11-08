/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: spinwalk.cu
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : Monte Carlo simulation of the spin dynamics in the presence of off-resonance fields.
 * -------------------------------------------------------------------------- */

// compile(lin) :  nvcc ./src/spinwalk.cu ./src/kernels.cu -I ./include/ -Xptxas -v -O3  -arch=compute_75 -code=sm_75  -Xcompiler -fopenmp -o spinwalk
// compile(win) :  nvcc ./src/spinwalk.cu ./src/kernels.cu -I ./include/ -Xptxas -v -O3  -arch=compute_86 -code=sm_86  -Xcompiler /openmp -std=c++17 -o spinwalk

#include <chrono>
#include <iomanip>
#include <algorithm>
#include <execution>
#include <filesystem>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include "indicators.hpp"
#include "CLI11.hpp"
#include "version.h"
#include "file_utils.h"
#include "shapes/cylinder.h"
#include "shapes/sphere.h"
#include "kernels.cuh"

#ifdef __CUDACC__
#include "helper.cuh"
#include <cuda_runtime.h>
#include "helper_cuda.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/iterator/constant_iterator.h>
#endif

#define BLOCKS  256
#define ERR_MSG  "\033[1;31mError:\033[0m "

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
    size_t numGrid = (param.n_spins + BLOCKS - 1) / BLOCKS;
    
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
            BOOST_LOG_TRIVIAL(info) << "GPU " << 1 << ") Generating random initial position for spins... (seed = " << param.seed << ", GPU grid = " << numGrid << " x " << BLOCKS << ")";
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
                        BOOST_LOG_TRIVIAL(warning) << "Virtual FoV (= " << param_local.fov[i] << ") is smaller than minimum convergence length (= " << min_convergence << ")!";
                        BOOST_LOG_TRIVIAL(warning) << "Original FoV (= " << param.fov[i] << ") FoV Scale (= " << fov_scale[sl]  << ")!";
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
            cu_sim<<<numGrid, BLOCKS, 0>>>(d_param, thrust::raw_pointer_cast(d_pFieldMap.data()), 
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
            bar.set_progress(100 * (sl+1)/float(param.n_fov_scale));
        }
        
        auto end_run     = std::chrono::high_resolution_clock::now();        
        auto elapsed_sim    = std::chrono::duration_cast<std::chrono::milliseconds>(end_run - start_sim).count() / 1000.0;
        auto elapsed_config = std::chrono::duration_cast<std::chrono::milliseconds>(end_run - start_config).count() / 1000.0;
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

        std::vector<float> TE;
        for(int i=0; i<param.n_TE; i++) TE.push_back(param.TE_us[i]*param.timestep_us*1e-6); 
        dims[0] = TE.size(); dims[1] = 1; dims[2] = 1; dims[3] = 1;
        file_utils::save_h5(f, TE.data(), dims, "TE", "float");
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
    bool bStatus = true, phantom_cylinder = false, phantom_sphere = false;
    std::string phantom_output;
    float phantom_radius=50.f, phantom_fov=1000.f, phantom_dchi=0.11e-6, phantom_oxy_level=0.75, phantom_orientation=90.f, phantom_volume_fraction=4.f;
    int32_t phantom_resolution=500, phantom_seed=-1, device_id = 0;
    std::vector<std::string>  config_files;
    std::vector<float> bvalue;
    
    // ========== parse command line arguments ==========
    CLI::App app{""};
    app.set_version_flag("-v,--version", get_verion());

    auto subcommand_phantom = app.add_subcommand("phantom", "Generate Numerical Phantom");
    subcommand_phantom->add_flag("-c,--cylinder", phantom_cylinder, "Fill phantom with cylinders");
    subcommand_phantom->add_flag("-s,--sphere", phantom_sphere, "Fill phantom with spheres");
    subcommand_phantom->add_option("-r,--radius", phantom_radius, "Radius of the cylinders/spheres in um (negative value = random radius)")->capture_default_str();
    subcommand_phantom->add_option("-n,--orientation", phantom_orientation, "Orientation of the cylinders in degree with respect to B0")->capture_default_str();
    subcommand_phantom->add_option("-v,--volume_fraction", phantom_volume_fraction, "Fraction of shapes volume to FoV volume <0.0 100.0>")->check(CLI::Range(0.0, 100.0))->capture_default_str();
    subcommand_phantom->add_option("-f,--fov", phantom_fov, "Voxel field of view in um (isotropic)")->capture_default_str()->check(CLI::PositiveNumber);
    subcommand_phantom->add_option("-z,--resolution", phantom_resolution, "Base resolution")->capture_default_str()->check(CLI::PositiveNumber);
    subcommand_phantom->add_option("-d,--dchi", phantom_dchi, "Susceptibility difference between fully deoxygenated blood and tissue (default: 0.11e-6 in cgs units)")->capture_default_str();
    subcommand_phantom->add_option("-y,--oxy_level", phantom_oxy_level, "Blood oxygenetation level <0.0 1.0> (-1 = exclude off-resonance effect and only generate the mask)")->capture_default_str();
    subcommand_phantom->add_option("-e,--seed", phantom_seed, "Seed for random number generator in phantom creator (-1 = random seed)")->capture_default_str();
    subcommand_phantom->add_option("-o,--output_file", phantom_output, "Path to save phantom (h5 format)")->capture_default_str();

    auto subcommand_sim = app.add_subcommand("sim", "Run Monte-Carlo Simulation");
    subcommand_sim->add_option("-c,--configs", config_files, "Config. files as many as you want. e.g. -c config1.ini config2.ini ... configN.ini")->check(CLI::ExistingFile);
    subcommand_sim->add_flag("-g,--gen_def_config", "Generate a default configuration file and store in the current folder");
#ifdef __CUDACC__
    subcommand_sim->add_flag("-p,--use_cpu", "Only run on CPU (default: GPU)");
    subcommand_sim->add_option("-d,--device", device_id, "Select GPU device (if there are multiple GPUs)");

    auto callback_gpu_info = [](int count){print_device_info();  exit(0);};
    app.add_flag("-g,--gpu_info", callback_gpu_info, "Print GPU information");
#endif

    auto subcommand_config = app.add_subcommand("config", "Generate Configuration File");
    subcommand_config->add_flag("-d,--default", "Generate a default configuration file and store in the current folder");
    subcommand_config->add_flag("-g,--gre", "Generate a configuration file for GRE sequence");
    subcommand_config->add_flag("-s,--se", "Generate a configuration file for Spin-Echo sequence");
    subcommand_config->add_flag("-b,--bssfp", "Generate a configuration file for bSSFP sequence");
    subcommand_config->add_option("-p,--pgse", bvalue, "Generate a configuration file for Pulsed-Gradient Spin-Echo sequence. The switch is followed by b-value, \xCE\xB4, and \xCE\x94. e.g. -p 1000 20 20")->expected(3);

    CLI11_PARSE(app, argc, argv);
    if(app.count_all() == 1)
    {
        std::cout << app.help() << '\n';
        return 0;
    }
    if (subcommand_phantom->parsed() && phantom_cylinder == phantom_sphere)
    {  
        if (phantom_cylinder == true)      
            std::cout << "Error! Please select either --cylinder or --sphere, not both!"<< '\n';
        std::cout << subcommand_phantom->help() << '\n';
        return 0;
    }
    if (subcommand_sim->parsed() && config_files.size() == 0 && subcommand_sim->count("--gen_def_config") == 0)
    {
        std::cout << subcommand_sim->help() << '\n';
        return 0;
    }
    if (subcommand_config->parsed())
    {
        std::cout << subcommand_config->help() << '\n';
        return 0;
    }

    // ========== logo ==========
    print_logo();    

    // ========== setup log ==========
    std::string log_filename = "spinwalk_" + std::to_string(device_id) + ".log";
    auto fileSink = bl::add_file_log(bl::keywords::file_name=log_filename, bl::keywords::target_file_name = log_filename, bl::keywords::format = "[%TimeStamp%] [%Severity%]: %Message%", bl::keywords::auto_flush = true);
    bl::add_common_attributes();
    std::cout << "Log file location: " << std::filesystem::current_path() / log_filename << '\n';

    // ========== generate phantom ==========
    if (phantom_cylinder)
    {
        cylinder cyl(phantom_fov, phantom_resolution, phantom_dchi, phantom_oxy_level, phantom_radius, phantom_volume_fraction, phantom_orientation, phantom_seed, phantom_output);
        cyl.run();
    }
    if (phantom_sphere)
    {
        sphere sph(phantom_fov, phantom_resolution, phantom_dchi, phantom_oxy_level, phantom_radius, phantom_volume_fraction, phantom_seed, phantom_output);
        sph.run();
    }

    // ========== generate default config ==========
    if (subcommand_sim->count("--gen_def_config"))
    {
        if(generate_default_config("./default_config.ini") == false)
            return 0;
        std::cout << "Default configuration file is generated in the current folder." << '\n';
    }

    // ========== loop over configs and simulate ==========
    if (config_files.size() == 0)
        return 0;
#ifdef __CUDACC__    
    if (subcommand_sim->count("--use_cpu") == 0)
        if(check_CUDA() == false)
            return 0;
#endif   

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
#ifdef __CUDACC__
        param.no_gpu = subcommand_sim->count("--use_cpu");
#else
        param.no_gpu = true;
#endif
        
        // ========== read config file ==========
        bStatus &= file_utils::read_config(cfile, &param, fov_scale, filenames);
            
        // ========== dump settings ==========
        bStatus &= dump_settings(param, filenames, fov_scale);

        // ========== Check GPU memory ==========
#ifdef __CUDACC__
        bStatus &= check_memory_size(param.get_required_memory(getDeviceCount()));
#endif
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
