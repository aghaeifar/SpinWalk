#ifndef SIM_H
#define SIM_H

#include <string>
#include "host_device_array.cuh"

struct simulation_parameters;

namespace sim {
class config_reader;

class monte_carlo {
    public:
        monte_carlo();
        virtual ~monte_carlo();
        virtual bool run(std::string config_filename);
        virtual void save(std::string filename);
#ifdef __CUDACC__
        void set_gpu_disabled(bool disabled){gpu_disabled = disabled;}
#endif    

    protected:
        virtual void allocate_memory();
        virtual bool read_phantom(std::string filename);
        virtual bool initialize_position(std::string filename, size_t seed);
        virtual bool initialize_magnetization(std::string filename);
        virtual size_t get_total_memory() const;

    private:        
        host_device_array<float>      fieldmap;        
        host_device_array<float>      XYZ0;      // memory layout(row-major): [n_spins x 3]
        host_device_array<float>      XYZ0_scaled;     // memory layout(row-major): [n_spins x 3]
        host_device_array<float>      XYZ1;      // memory layout(row-major): [n_fov_scale x n_spins x timepoints x 3] or [n_fov_scale x n_spins x 1 x 3]
        host_device_array<float>      M0;        // memory layout(row-major): [n_spins x 3]
        host_device_array<float>      M1;        // memory layout(row-major): [n_fov_scale x n_spins x n_TE x 3]
        host_device_array<uint8_t>    T;         // memory layout(row-major): [n_fov_scale x n_spins x n_TE x 1]
        host_device_array<uint8_t>    mask;
        host_device_array<float>      fov;
        simulation_parameters         *param;
        config_reader                 *config;
#ifdef __CUDACC__
        bool gpu_disabled;
        int32_t device_count;
#endif

};



} // namespace sim


#endif // SIM_H