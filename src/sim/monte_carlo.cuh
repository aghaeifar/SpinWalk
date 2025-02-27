#ifndef SIM_H
#define SIM_H

#include <string>
#include <vector>
#include "config_reader.h"
#include "simulation_parameters.cuh"

namespace sim {
class monte_carlo {
    public:
        monte_carlo(bool gpu_disabled=false, int32_t device_id=0);
        virtual ~monte_carlo();
        virtual bool run(std::string config_filename);
        virtual void save(std::string filename);

    protected:
        virtual void allocate_memory();
        virtual bool read_phantom(std::string filename);
        virtual bool initialize_position(std::string filename, size_t seed);
        virtual bool initialize_magnetization(std::string filename);
        virtual size_t get_total_memory() const;
        #ifdef __CUDACC__
        virtual void cleanup_device();
        #endif


    private:        
        std::vector<float>      fieldmap;        
        std::vector<float>      XYZ0;      // memory layout(row-major): [n_spins x 3]
        std::vector<float>      XYZ0_scaled;     // memory layout(row-major): [n_spins x 3]
        std::vector<float>      XYZ1;      // memory layout(row-major): [n_fov_scale x n_spins x timepoints x 3] or [n_fov_scale x n_spins x 1 x 3]
        std::vector<float>      M0;        // memory layout(row-major): [n_spins x 3]
        std::vector<float>      M1;        // memory layout(row-major): [n_fov_scale x n_spins x n_TE x 3]
        std::vector<uint8_t>    T;         // memory layout(row-major): [n_fov_scale x n_spins x n_TE x 1]
        std::vector<uint8_t>    mask;
        std::vector<float>      fov;
        config_reader           config;
        parameters              param;
        parameters_hvec         param_hvec;
        parameters_uvec         param_uvec;
        bool gpu_disabled;
        
#ifdef __CUDACC__
        parameters_dvec         param_dvec;
        thrust::device_vector<float>      d_fieldmap;        
        thrust::device_vector<float>      d_XYZ0;      // memory layout(row-major): [n_spins x 3]
        thrust::device_vector<float>      d_XYZ1;      // memory layout(row-major): [n_fov_scale x n_spins x timepoints x 3] or [n_fov_scale x n_spins x 1 x 3]
        thrust::device_vector<float>      d_M0;        // memory layout(row-major): [n_spins x 3]
        thrust::device_vector<float>      d_M1;        // memory layout(row-major): [n_fov_scale x n_spins x n_TE x 3]
        thrust::device_vector<uint8_t>    d_T;         // memory layout(row-major): [n_fov_scale x n_spins x n_TE x 1]
        thrust::device_vector<uint8_t>    d_mask;
#endif

};



} // namespace sim


#endif // SIM_H