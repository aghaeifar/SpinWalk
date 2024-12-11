// standard libraries
#include <chrono>
#include <execution>
#include <algorithm>
#include <filesystem>

// custom headers
#include "monte_carlo.cuh"
#include "kernels.cuh"
#include "h5_helper.h"
#include "barkeep.h"
#include "config_reader.h"
#include "simulation_parameters.cuh"
#include "definitions.h"

// boost headers
#include <boost/log/trivial.hpp>

// CUDA libraries
#ifdef __CUDACC__
#include "device_helper.cuh"
#include <cuda_runtime.h>
#endif


#define BLOCKS  256

namespace bl = boost::log;

namespace sim
{

monte_carlo::monte_carlo() 
{
    param  = new simulation_parameters();
    config = new config_reader();
#ifdef __CUDACC__
    gpu_disabled  = false;
    checkCudaErrors(cudaGetDeviceCount(&device_count));
    BOOST_LOG_TRIVIAL(info) << "Number of available GPU(s): " << device_count; 
#endif    
}

monte_carlo::~monte_carlo()
{
    if(param)  delete param;
    if(config) delete config;
}

void monte_carlo::allocate_memory()
{
    size_t trj_dim_size  = param->enRecordTrajectory ? param->n_timepoints * (param->n_dummy_scan + 1) : 1;
    XYZ0.data.resize(param->n_spins * 3);     // memory layout(row-major): [n_spins x 3]
    XYZ0_scaled.data.resize(XYZ0.data.size());       // memory layout(row-major): [n_spins x 3]
    XYZ1.data.resize(trj_dim_size * param->n_scales * XYZ0.data.size());     // memory layout(row-major): [n_scale x n_spins x timepoints x 3] or [n_scale x n_spins x 1 x 3]
    M0.data.resize(param->n_spins * 3);       // memory layout(row-major): [n_spins x 3]
    M1.data.resize(param->n_scales * param->TE_us.data.size() * M0.data.size());    // memory layout(row-major): [n_scale x n_spins x n_TE x 3]
    T.data.resize(M1.data.size()/3);                 // memory layout(row-major): [n_scale x n_spins x n_TE x 1]
    gradient_mTm.data.resize(param->gradient_mTm.data.size());
}

size_t monte_carlo::get_total_memory() const
{
    size_t total_memory = 0;
    total_memory += XYZ0_scaled.data.size()  * sizeof(float);
    total_memory += XYZ1.data.size()         * sizeof(float);
    total_memory += M0.data.size()           * sizeof(float);
    total_memory += M1.data.size()           * sizeof(float);
    total_memory += T.data.size()            * sizeof(uint8_t);
    total_memory += mask.data.size()         * sizeof(uint8_t);
    total_memory += fieldmap.data.size()     * sizeof(float);
    return total_memory >> 20; // convert to MB
}

bool monte_carlo::read_phantom(std::string filename)
{
    fov.data.resize(3);
    if(h5_helper::read(filename, "fieldmap", true, fieldmap.data) == false)
        fieldmap.data.clear();
    if(h5_helper::read(filename, "mask", true, mask.data) == false)
        return false;   
    if(h5_helper::read(filename, "fov", false, fov.data) == false)
        return false;   

    std::vector<size_t> phantom_size;  
    if(h5_helper::size(filename, "mask", phantom_size) == false)
        return false;
    std::copy(phantom_size.begin(), phantom_size.end(), param->phantom_size);

    uint32_t n_substrate = *std::max_element(std::execution::par, mask.data.begin(), mask.data.end()) + 1;
    if (n_substrate > param->n_substrate)
    {
        BOOST_LOG_TRIVIAL(error) << "The number of substrate types in the mask does not match the number of substrate types in the config file: " << n_substrate << " vs " << param->n_substrate;
        return false;
    }

    BOOST_LOG_TRIVIAL(info) << "Size = " << phantom_size[0] << " x " << phantom_size[1] << " x " << phantom_size[2] << std::endl;
    BOOST_LOG_TRIVIAL(info) << "FoV = " << fov.data[0]*1e6 << " x " << fov.data[1]*1e6 << " x " << fov.data[2]*1e6 << " um^3" << std::endl;
    return true;
}

bool monte_carlo::initialize_position(std::string filename, size_t seed)
{
    BOOST_LOG_TRIVIAL(info) << "Initializing positions...";
    if (filename.empty() == false){
        BOOST_LOG_TRIVIAL(info) << "Reading initial positions from file: " << filename;
        if(h5_helper::read(filename, "XYZ", false, XYZ0.data) == false)
            return false; 
        // check values from file are in FoV
        uint32_t ind=0;
        if (std::any_of(std::execution::par, XYZ0.data.begin(), XYZ0.data.end(), [this, &ind](float x){return x < 0 || x > fov.data[ind++%3];})) {
            BOOST_LOG_TRIVIAL(error) << "Initial positions are outside the FoV.";
            return false;
        }
        return true;
    }
    // if no filename is provided, generate random positions
    BOOST_LOG_TRIVIAL(info) << "Generating random positions within 98% of the FoV with seed = " << seed;
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dist_initial_x(0.01*fov.data[0], 0.99*fov.data[0]);
    std::uniform_real_distribution<float> dist_initial_y(0.01*fov.data[1], 0.99*fov.data[1]);
    std::uniform_real_distribution<float> dist_initial_z(0.01*fov.data[2], 0.99*fov.data[2]);

    for (size_t i = 0; i < XYZ0.data.size() / 3; i++){
        XYZ0.data[3*i+0] = dist_initial_x(gen);
        XYZ0.data[3*i+1] = dist_initial_y(gen);
        XYZ0.data[3*i+2] = dist_initial_z(gen);
    }    
    return true;
}

bool monte_carlo::initialize_magnetization(std::string filename)
{
    BOOST_LOG_TRIVIAL(info) << "Initializing magnetization...";
    if (filename.empty() == false)
        if(h5_helper::read(filename, "M", false, M0.data) == false)
            return false; 

    BOOST_LOG_TRIVIAL(info) << "Generating M0(0, 0, 1)..." << std::endl;
    uint32_t index = 0;
    std::generate(M0.data.begin(), M0.data.end(), [&index](){return (index++ % 3 == 2) ? 1.f : 0.f;});
    return true;
}

void monte_carlo::save(std::string filename)
{
    M1.copy_to_host();
    std::vector<size_t> dims = {param->n_scales, param->n_spins, param->TE_us.data.size(), 3};
    h5_helper::write(filename, "M", dims, M1.data);

    XYZ1.copy_to_host();
    dims[2] = param->enRecordTrajectory ? param->n_timepoints * (param->n_dummy_scan + 1) : 1;
    h5_helper::write(filename, "XYZ", dims, XYZ1.data);

    T.copy_to_host();
    dims[3] = 1; dims[2] = param->TE_us.data.size();
    h5_helper::write(filename, "T", dims, T.data);

    dims[0] = config->get_scales().size(); dims[1] = 1; dims[2] = 1; dims[3] = 1;
    h5_helper::write(filename, "scales", dims, config->get_scales());

    std::vector<float> TE_us;
    for(int i=0; i<param->TE_us.data.size(); i++) TE_us.push_back(param->TE_us.data[i]*param->timestep_us*1e-6); 
    dims[0] = TE_us.size(); dims[1] = 1; dims[2] = 1; dims[3] = 1;
    h5_helper::write(filename, "TE", dims, TE_us);
}

bool monte_carlo::run(std::string config_filename) // simulation_parameters param, std::map<std::string, std::vector<std::string> > filenames, std::vector<double> scale
{   
    auto start_run = std::chrono::high_resolution_clock::now();
    // ========== read config file ==========
    if(config->prepare(config_filename, param) == false)
        return false; 
    if (param->prepare() == false)
        return false;
    allocate_memory();
    
    BOOST_LOG_TRIVIAL(info) << "\n" << std::string(20, '-') << "\nSimulation parameters:\n" << param->dump() << "\n" << std::string(20, '-') ;
    size_t trj  = param->enRecordTrajectory ? param->n_timepoints * (param->n_dummy_scan + 1) : 1;
    size_t ind_fieldmap = 0;
    gradient_mTm.data = param->gradient_mTm.data;
    for (auto &file_phantom : config->get_filename("PHANTOM")){
        BOOST_LOG_TRIVIAL(info) << "Simulating phantom: " << file_phantom;
        if(read_phantom(file_phantom) == false)
            return false;        
        if(initialize_position(config->get_filename("XYZ0")[ind_fieldmap], param->seed) == false)
            return false;
        if(initialize_magnetization(config->get_filename("M0")[ind_fieldmap]) == false)
            return false;

        XYZ0_scaled.data      = XYZ0.data;
        param->matrix_length  = mask.data.size(); // update the matrix length based on the mask size from the recent read
        param->fieldmap_exist = fieldmap.data.size() > 0;

        // convert fieldmap from T to degree per timestep
        float Tesla2deg_pertimestep = param->B0 * param->timestep_us * 1e-6 * GAMMA * RAD2DEG;
        BOOST_LOG_TRIVIAL(info) << "Conversion factor from T to degree per timestep: " << Tesla2deg_pertimestep;
        if(param->fieldmap_exist) 
            std::transform(std::execution::par_unseq, fieldmap.data.begin(), fieldmap.data.end(), fieldmap.data.begin(), [Tesla2deg_pertimestep](auto x) { return x*Tesla2deg_pertimestep;});

        // ========== move to GPU memory ========== 
        bool device = CPU;
#ifdef __CUDACC__
        if (gpu_disabled == false) { 
            BOOST_LOG_TRIVIAL(info) << "Moving data to GPU memory.";
            // calculate required memory and avaialbe memory 
            if (check_memory_size(get_total_memory()) == false)
                return false; 
            device = GPU;
        }
#endif

        XYZ1.init(device, false);
        M0.init(device, true); 
        M1.init(device, false);
        T.init(device, false); 
        XYZ0_scaled.init(device, true);
        XYZ0.init(device, true); 
        mask.init(device, true); 
        fieldmap.init(device, true);
        gradient_mTm.init(device, false);
        param->init(device, true); 
        for (int i = 0; i < 3; i++) // FoV scaling
            param->fov[i] = fov.data[i]; 
        // ========== run ==========   
        uint32_t ind_scale = 0;
        std::vector<uint32_t> v(param->n_spins);
        std::generate(std::execution::seq, v.begin(), v.end(), [n = 0] () mutable { return n++; }); 
        
        auto start_sim = std::chrono::high_resolution_clock::now();
        auto bar = barkeep::ProgressBar(&ind_scale, {.total = param->n_scales, .message = "Simulating", .style = barkeep::ProgressBarStyle::Rich,});
        for (const auto scale : config->get_scales())
        {   
            BOOST_LOG_TRIVIAL(info) << "Simulating scale " << scale; 
            if(config->get_scale_type() == e_scale_type::s_fov){           
                std::transform(std::execution::par_unseq, XYZ0.data.begin(), XYZ0.data.end(), XYZ0_scaled.data.begin(), [scale](auto& c){return c*scale;}); 
                for (int i = 0; i < 3; i++) // FoV scaling
                    param->fov[i] = scale * fov.data[i]; 
                XYZ0_scaled.init(device, true);
             } else if (config->get_scale_type() == e_scale_type::s_gradient){                                
                std::transform(std::execution::par_unseq, gradient_mTm.data.begin(), gradient_mTm.data.end(), param->gradient_mTm.data.begin(), [scale](auto& c){return c*scale;}); 
                param->gradient_mTm.init(device, true);
             }

            // here we need to check voxel size and step size to make sure that the simulation is stable: doi:10.1016/j.neuroimage.2018.06.046 & https://submissions.mirasmart.com/ISMRM2024/Itinerary/PresentationDetail.aspx?evdid=4684
            
            // ========== simulation kernel  ==========
#ifdef __CUDACC__
            if(gpu_disabled){
#endif           
                std::for_each(std::execution::par_unseq, v.begin(), v.end(), [&](int spin) {sim(param, fieldmap.ptr, mask.ptr, M0.ptr, XYZ0_scaled.ptr, 
                                                                                                M1.ptr + 3*param->TE_us.data.size()*param->n_spins*ind_scale, 
                                                                                                XYZ1.ptr + 3*param->n_spins*trj*ind_scale, 
                                                                                                T.ptr + param->TE_us.data.size()*param->n_spins*ind_scale,
                                                                                                spin);});
#ifdef __CUDACC__  
            }else{         
                size_t numGrid = (param->n_spins + BLOCKS - 1) / BLOCKS;
                cu_sim<<<numGrid, BLOCKS, 0>>>(param, fieldmap.ptr, mask.ptr, M0.ptr, XYZ0_scaled.ptr, 
                                               M1.ptr + 3*param->TE_us.data.size()*param->n_spins*ind_scale, 
                                               XYZ1.ptr + 3*param->n_spins*trj*ind_scale, 
                                               T.ptr + param->TE_us.data.size()*param->n_spins*ind_scale); 
                gpuCheckKernelExecutionError(__FILE__, __LINE__);
            }
#endif     
            ind_scale++;       
        }        
        bar->done();

        auto elapsed_sim = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_sim).count() / 1000.0;
        int precision = elapsed_sim>10 ? 0 : (elapsed_sim > 1 ? 1 : 3);
        BOOST_LOG_TRIVIAL(info) << "Simulation took " << std::fixed << std::setprecision(precision) <<  elapsed_sim << " seconds.";

        // ========== save results ========== 
        BOOST_LOG_TRIVIAL(info) << "Saving the results to disk.";
        save(config->get_output_filename(ind_fieldmap));
        ind_fieldmap++;
    }

    auto elapsed_run = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_run).count() / 1000.0;
    int precision = elapsed_run>10 ? 0 : (elapsed_run > 1 ? 1 : 3);
    BOOST_LOG_TRIVIAL(info) << "Entire run took " << std::fixed << std::setprecision(precision) <<  elapsed_run << " seconds.";
    return true;
} 

} // namespace sim