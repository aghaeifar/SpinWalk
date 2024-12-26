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

#include "definitions.h"

// boost headers
#include <boost/log/trivial.hpp>

// CUDA libraries
#ifdef __CUDACC__
#include "device_helper.cuh"
#include <cuda_runtime.h>
#include "helper_cuda.h"
#endif


#define BLOCKS  256

namespace bl = boost::log;

namespace sim
{

monte_carlo::monte_carlo(bool gpu_disabled, int32_t device_id)
{
    this->gpu_disabled = true;
#ifdef __CUDACC__
    this->gpu_disabled  = gpu_disabled;
    if(gpu_disabled == false){
        if(sim::check_CUDA() == false){
            std::cout << WARN_MSG << "No GPU Device found! switching to CPU mode." << std::endl;
            this->gpu_disabled = gpu_disabled = true;
        }
    }
    if(gpu_disabled == false){
        uint32_t device_count = sim::get_device_count();
        if (device_id >= device_count){
            std::cout << ERR_MSG << "Device ID " << device_id << " is not available! Number of available GPU(s) is " << device_count << " ,switching to CPU mode!" << std::endl;
            this->gpu_disabled = gpu_disabled = true;
        } else {
            BOOST_LOG_TRIVIAL(info) << "Number of available GPU(s): " << device_count; 
            cudaSetDevice(device_id);
        }
    }
#endif  
}

monte_carlo::~monte_carlo()
{
}

void monte_carlo::allocate_memory()
{
    size_t trj_dim_size  = param.enRecordTrajectory ? param.n_timepoints * (param.n_dummy_scan + 1) : 1;
    XYZ0.resize(param.n_spins * 3);     // memory layout(row-major): [n_spins x 3]
    XYZ0_scaled.resize(XYZ0.size());       // memory layout(row-major): [n_spins x 3]
    XYZ1.resize(trj_dim_size * param.n_scales * XYZ0.size());     // memory layout(row-major): [n_scale x n_spins x timepoints x 3] or [n_scale x n_spins x 1 x 3]
    M0.resize(param.n_spins * 3);       // memory layout(row-major): [n_spins x 3]
    M1.resize(param.n_scales * param_hvec.TE_us.size() * M0.size());    // memory layout(row-major): [n_scale x n_spins x n_TE x 3]
    T.resize(M1.size()/3);                 // memory layout(row-major): [n_scale x n_spins x n_TE x 1]
}

size_t monte_carlo::get_total_memory() const
{
    size_t total_memory = 0;
    total_memory += XYZ0_scaled.size()  * sizeof(float);
    total_memory += XYZ1.size()         * sizeof(float);
    total_memory += M0.size()           * sizeof(float);
    total_memory += M1.size()           * sizeof(float);
    total_memory += T.size()            * sizeof(uint8_t);
    total_memory += mask.size()         * sizeof(uint8_t);
    total_memory += fieldmap.size()     * sizeof(float);
    return total_memory >> 20; // convert to MB
}

bool monte_carlo::read_phantom(std::string filename)
{
    fov.resize(3);
    if(h5_helper::read(filename, "fieldmap", true, fieldmap) == false)
        fieldmap.clear();
    if(h5_helper::read(filename, "mask", true, mask) == false)
        return false;   
    if(h5_helper::read(filename, "fov", false, fov) == false)
        return false;   

    std::vector<size_t> phantom_size;  
    if(h5_helper::size(filename, "mask", phantom_size) == false)
        return false;
    std::copy(phantom_size.begin(), phantom_size.end(), param.phantom_size);

    uint32_t n_substrate = *std::max_element(std::execution::par, mask.begin(), mask.end()) + 1;
    if (n_substrate > param.n_substrate)
    {
        BOOST_LOG_TRIVIAL(error) << "The number of substrate types in the mask does not match the number of substrate types in the config file: " << n_substrate << " vs " << param.n_substrate;
        return false;
    }

    BOOST_LOG_TRIVIAL(info) << "Size = " << phantom_size[0] << " x " << phantom_size[1] << " x " << phantom_size[2] << std::endl;
    BOOST_LOG_TRIVIAL(info) << "FoV = " << fov[0]*1e6 << " x " << fov[1]*1e6 << " x " << fov[2]*1e6 << " um^3" << std::endl;
    return true;
}

bool monte_carlo::initialize_position(std::string filename, size_t seed)
{
    BOOST_LOG_TRIVIAL(info) << "Initializing positions...";
    if (filename.empty() == false){
        BOOST_LOG_TRIVIAL(info) << "Reading initial positions from file: " << filename;
        if(h5_helper::read(filename, "XYZ", false, XYZ0) == false)
            return false; 
        // check values from file are in FoV
        uint32_t ind=0;
        if (std::any_of(std::execution::par, XYZ0.begin(), XYZ0.end(), [this, &ind](float x){return x < 0 || x > fov[ind++%3];})) {
            BOOST_LOG_TRIVIAL(error) << "Initial positions are outside the FoV.";
            return false;
        }
        return true;
    }
    // if no filename is provided, generate random positions
    BOOST_LOG_TRIVIAL(info) << "Generating random positions within 98% of the FoV with seed = " << seed;
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dist_initial_x(0.01*fov[0], 0.99*fov[0]);
    std::uniform_real_distribution<float> dist_initial_y(0.01*fov[1], 0.99*fov[1]);
    std::uniform_real_distribution<float> dist_initial_z(0.01*fov[2], 0.99*fov[2]);

    for (size_t i = 0; i < XYZ0.size() / 3; i++){
        XYZ0[3*i+0] = dist_initial_x(gen);
        XYZ0[3*i+1] = dist_initial_y(gen);
        XYZ0[3*i+2] = dist_initial_z(gen);
    }    
    return true;
}

bool monte_carlo::initialize_magnetization(std::string filename)
{
    BOOST_LOG_TRIVIAL(info) << "Initializing magnetization...";
    if (filename.empty() == false)
        if(h5_helper::read(filename, "M", false, M0) == false)
            return false; 

    BOOST_LOG_TRIVIAL(info) << "Generating M0(0, 0, 1)..." << std::endl;
    uint32_t index = 0;
    std::generate(M0.begin(), M0.end(), [&index](){return (index++ % 3 == 2) ? 1.f : 0.f;});
    return true;
}

void monte_carlo::save(std::string filename)
{
#ifdef __CUDACC__
    if(gpu_disabled == false){
        thrust::copy(d_M1.begin(), d_M1.end(), M1.begin());
        thrust::copy(d_XYZ1.begin(), d_XYZ1.end(), XYZ1.begin());
        thrust::copy(d_T.begin(), d_T.end(), T.begin());
    }
#endif

    std::vector<size_t> dims = {param.n_scales, param.n_spins, param_hvec.TE_us.size(), 3};
    h5_helper::write(filename, "M", dims, M1);

    dims[2] = param.enRecordTrajectory ? param.n_timepoints * (param.n_dummy_scan + 1) : 1;
    h5_helper::write(filename, "XYZ", dims, XYZ1);

    dims[3] = 1; dims[2] = param_hvec.TE_us.size();
    h5_helper::write(filename, "T", dims, T);

    dims[0] = config.get_scales().size(); dims[1] = 1; dims[2] = 1; dims[3] = 1;
    h5_helper::write(filename, "scales", dims, config.get_scales());

    std::vector<float> TE_us;
    for(int i=0; i<param_hvec.TE_us.size(); i++) TE_us.push_back(param_hvec.TE_us[i]*param.timestep_us*1e-6); 
    dims[0] = TE_us.size(); dims[1] = 1; dims[2] = 1; dims[3] = 1;
    h5_helper::write(filename, "TE", dims, TE_us);
}

bool monte_carlo::run(std::string config_filename) // simulation_parameters param, std::map<std::string, std::vector<std::string> > filenames, std::vector<double> scale
{   
    auto start_run = std::chrono::high_resolution_clock::now();
    // ========== read config file ==========
    if(config.prepare(config_filename, &param, &param_hvec) == false)
        return false; 
    if (param.prepare(param_hvec) == false)
        return false;
    allocate_memory();
    
    // BOOST_LOG_TRIVIAL(info) << "\n" << std::string(20, '-') << "\nSimulation parameters:\n" << param.dump() << "\n" << std::string(20, '-') ;
    size_t trj  = param.enRecordTrajectory ? param.n_timepoints * (param.n_dummy_scan + 1) : 1;
    size_t ind_fieldmap = 0;
    std::vector<float> gradient_mTm_orig = param_hvec.gradient_mTm;

    param_uvec.copy_from_host(param_hvec);
#ifdef __CUDACC__
    if (gpu_disabled == false) { 
        param_dvec.copy_from_host(param_hvec);
        param_uvec.copy_from_device(param_dvec);
    }
#endif
    for (auto &file_phantom : config.get_filename("PHANTOM")){
        BOOST_LOG_TRIVIAL(info) << "Simulating phantom: " << file_phantom;
        if(read_phantom(file_phantom) == false)
            return false;        
        if(initialize_position(config.get_filename("XYZ0")[ind_fieldmap], param.seed) == false)
            return false;
        if(initialize_magnetization(config.get_filename("M0")[ind_fieldmap]) == false)
            return false;

        XYZ0_scaled      = XYZ0;
        param.matrix_length  = mask.size(); // update the matrix length based on the mask size from the recent read
        param.fieldmap_exist = fieldmap.size() > 0;
       
        // convert fieldmap from T to degree per timestep
        float Tesla2deg_pertimestep = param.B0 * param.timestep_us * 1e-6 * GAMMA * RAD2DEG;
        BOOST_LOG_TRIVIAL(info) << "Conversion factor from T to degree per timestep: " << Tesla2deg_pertimestep;
        if(param.fieldmap_exist) 
            std::transform(std::execution::par_unseq, fieldmap.begin(), fieldmap.end(), fieldmap.begin(), [Tesla2deg_pertimestep](auto x) { return x*Tesla2deg_pertimestep;});

        // ========== move to GPU memory ========== 
#ifdef __CUDACC__
        if (gpu_disabled == false) { 
            BOOST_LOG_TRIVIAL(info) << "Moving data to GPU memory.";
            // calculate required memory and avaialbe memory 
            if (check_memory_size(get_total_memory()) == false)
                return false; 
            d_fieldmap = fieldmap;
            d_XYZ0 = XYZ0;
            d_XYZ1 = XYZ1;
            d_mask = mask;
            d_M0 = M0;
            d_M1 = M1;
            d_T = T;
        }
#endif

        for (int i = 0; i < 3; i++) // FoV scaling
            param.fov[i] = fov[i]; 
        // ========== run ==========   
        uint32_t ind_scale = 0;
        std::vector<uint32_t> v(param.n_spins);
        std::generate(std::execution::seq, v.begin(), v.end(), [n = 0] () mutable { return n++; }); 
        
        auto start_sim = std::chrono::high_resolution_clock::now();
        auto bar = barkeep::ProgressBar(&ind_scale, {.total = param.n_scales, .message = "Simulating", .style = barkeep::ProgressBarStyle::Rich,});
        for (const auto scale : config.get_scales())
        {   
            BOOST_LOG_TRIVIAL(info) << "Simulating scale " << scale; 
            // FoV scaling
            if(config.get_scale_type() == e_scale_type::s_fov){           
                std::transform(std::execution::par_unseq, XYZ0.begin(), XYZ0.end(), XYZ0_scaled.begin(), [scale](auto& c){return c*scale;}); 
                for (int i = 0; i < 3; i++) // FoV scaling
                    param.fov[i] = scale * fov[i]; 
#ifdef __CUDACC__
                if (gpu_disabled == false)
                    d_XYZ0 = XYZ0_scaled;
#endif
             } 
             // Gradient scaling
             else if (config.get_scale_type() == e_scale_type::s_gradient){                                
                std::transform(std::execution::par_unseq, gradient_mTm_orig.begin(), gradient_mTm_orig.end(), param_hvec.gradient_mTm.begin(), [scale](auto& c){return c*scale;}); 
#ifdef __CUDACC__
                if (gpu_disabled == false){
                    param_dvec.gradient_mTm = param_hvec.gradient_mTm;
                    param_uvec.gradient_mTm.ptr = thrust::raw_pointer_cast(param_dvec.gradient_mTm.data());
                }
#endif
             }

            // here we need to check voxel size and step size to make sure that the simulation is stable: doi:10.1016/j.neuroimage.2018.06.046 & https://submissions.mirasmart.com/ISMRM2024/Itinerary/PresentationDetail.aspx?evdid=4684
            
            // ========== simulation kernel  ==========
#ifdef __CUDACC__
            if(gpu_disabled){
#endif           
                std::for_each(std::execution::par_unseq, v.begin(), v.end(), [&](int spin) {sim(param, param_uvec, 
                                                                                            fieldmap.data(), 
                                                                                            mask.data(), 
                                                                                            M0.data(), 
                                                                                            XYZ0_scaled.data(), 
                                                                                            M1.data() + 3*param_hvec.TE_us.size()*param.n_spins*ind_scale, 
                                                                                            XYZ1.data() + 3*param.n_spins*trj*ind_scale, 
                                                                                            T.data() + param_hvec.TE_us.size()*param.n_spins*ind_scale,
                                                                                            spin);});
#ifdef __CUDACC__  
            }else{         
                size_t numGrid = (param.n_spins + BLOCKS - 1) / BLOCKS;
                cu_sim<<<numGrid, BLOCKS, 0>>>(param, param_uvec, 
                                                thrust::raw_pointer_cast(d_fieldmap.data()), 
                                                thrust::raw_pointer_cast(d_mask.data()),
                                                thrust::raw_pointer_cast(d_M0.data()),
                                                thrust::raw_pointer_cast(d_XYZ0.data()),
                                                thrust::raw_pointer_cast(d_M1.data() + 3*param_dvec.TE_us.size()*param.n_spins*ind_scale),
                                                thrust::raw_pointer_cast(d_XYZ1.data() + 3*param.n_spins*trj*ind_scale),
                                                thrust::raw_pointer_cast(d_T.data() + param_dvec.TE_us.size()*param.n_spins*ind_scale));
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
        save(config.get_output_filename(ind_fieldmap));
        ind_fieldmap++;
    }

    auto elapsed_run = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_run).count() / 1000.0;
    int precision = elapsed_run>10 ? 0 : (elapsed_run > 1 ? 1 : 3);
    BOOST_LOG_TRIVIAL(info) << "Entire run took " << std::fixed << std::setprecision(precision) <<  elapsed_run << " seconds.";
    return true;
} 

} // namespace sim