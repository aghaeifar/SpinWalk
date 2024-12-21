
#ifndef DEVICE_HELPER_CUH
#define DEVICE_HELPER_CUH

#include <iostream>
#include <boost/log/trivial.hpp>
#include "device_helper.cuh"

#ifdef __CUDACC__
#include "helper_cuda.h"
#endif

//---------------------------------------------------------------------------------------------
//  check for CUDA and GPU device
//---------------------------------------------------------------------------------------------

#define NO_CUDA_MSG "CUDA is not available in this build. Please recompile with CUDA enabled."

namespace sim
{
bool check_CUDA()
{
#ifdef __CUDACC__
    int32_t device = 0;
    cudaError_t error = cudaGetDeviceCount(&device);
    if (error != cudaSuccess)
    {
        std::cout << "\033[1;31mError:\033[0m " <<cudaGetErrorString(error) << "\n";
        return false;
    }
#endif
    return true;
}

uint32_t get_device_count()
{
    int32_t device_count = 0;
#ifdef __CUDACC__
    checkCudaErrors(cudaGetDeviceCount(&device_count));
#endif
    return device_count;
}

std::string convert_cuda_version(int32_t value) {
    int v = value/1000;
    int h = (value % 100) / 10;
    return std::to_string(v) + "." + std::to_string(h);
}

void print_device_info()
{
#ifdef __CUDACC__
    if (check_CUDA() == false)
        return;

    size_t free, total;    
    int32_t cuda_version, driver_version;  
    
    checkCudaErrors(cudaRuntimeGetVersion(&cuda_version));
    checkCudaErrors(cudaDriverGetVersion(&driver_version));
    std::cout << "The latest version of CUDA supported by the driver: "<< convert_cuda_version(driver_version) << ", current CUDA version: "<< convert_cuda_version(cuda_version) << "\n";

    int32_t device_count = get_device_count();
    std::cout <<"Number of devices: " << device_count << "\n";

    cudaDeviceProp device_properties;

    int device_id = 0;
    checkCudaErrors(cudaGetDevice(&device_id));
    cudaGetDeviceProperties(&device_properties, device_id);
    std::cout << device_properties.name << std::endl;
    std::cout << "-Compute Capability: " << device_properties.major << "."<< device_properties.minor << std::endl;
    // cudaDeviceReset();
    cudaMemGetInfo(&free, &total);
    std::cout << "-Free GPU Memory: " << (free >> 20) << " MB (out of " << (total >> 20) << " MB)" << std::endl;
#else
    std::cout << NO_CUDA_MSG << std::endl;
#endif
}
 
bool check_memory_size(size_t required_size_MB)
{
    bool memory_ok = true;
#ifdef __CUDACC__
    if (check_CUDA() == false)
        return false;

    size_t free, total;
    int32_t device_count = get_device_count();
    cudaDeviceProp device_properties;

    int device_id = 0;
    checkCudaErrors(cudaGetDevice(&device_id));
    checkCudaErrors(cudaGetDeviceProperties(&device_properties, device_id));
    checkCudaErrors(cudaMemGetInfo(&free, &total));
    free  = free  >> 20;
    total = total >> 20;
    BOOST_LOG_TRIVIAL(info) << device_properties.name << " " << device_id << ": total memroy= " << total << " MB, free memory= " << free << " MB, required memory= " << required_size_MB << " MB" << '\n';
    memory_ok = free<required_size_MB ? false:true;
    if(memory_ok == false)
        BOOST_LOG_TRIVIAL(fatal) << "Not enough GPU memory for the simulation in device! Required=" << required_size_MB <<" MB, Available=" << free << " MB";
#endif
    return memory_ok;
}

} // namespace sim
#endif // DEVICE_HELPER_CUH
