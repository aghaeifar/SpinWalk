
#ifndef _HELPER_H_
#define _HELPER_H_


#include <boost/log/trivial.hpp>
#include "helper_cuda.h"

//---------------------------------------------------------------------------------------------
//  check for CUDA and GPU device
//---------------------------------------------------------------------------------------------

bool check_CUDA()
{
    int32_t device = 0;
    cudaError_t error = cudaGetDeviceCount(&device);
    if (error != cudaSuccess)
    {
        std::cout << "\033[1;31mError:\033[0m " <<cudaGetErrorString(error) << "\n";
        return false;
    }
    return true;
}

uint32_t getDeviceCount()
{
    int32_t device_count = 0;
    checkCudaErrors(cudaGetDeviceCount(&device_count));
    return device_count;
}

void print_device_info()
{
    if (check_CUDA() == false)
        return;

    const int kb = 1024;
    const int mb = kb * kb;
    size_t free, total;    
    int32_t cuda_version, driver_version;  
    
    cudaRuntimeGetVersion(&cuda_version);
    cudaDriverGetVersion(&driver_version);
    std::cout << "Driver version: "<< driver_version << ", CUDA version: "<< cuda_version << "\n";

    int32_t device_count = getDeviceCount();
    std::cout <<"Number of devices: " << device_count << "\n";

    cudaDeviceProp device_properties;
    for(int i=0; i<device_count; i++)
    {
        cudaSetDevice(i);
        cudaGetDeviceProperties(&device_properties, i);
        std::cout << "Device " << i+1 << ", " << device_properties.name << std::endl;
        std::cout << "-Compute Capability: " << device_properties.major << "."<< device_properties.minor << std::endl;
        // cudaDeviceReset();
        cudaMemGetInfo(&free, &total);
        std::cout << "-Free GPU Memory: " << free / mb << " MB (out of " << total / mb << " MB)" << std::endl;
    }
}
 
bool check_memory_size(size_t required_size_MB)
{
    if (check_CUDA() == false)
        return false;

    size_t free, total;
    bool memory_ok = true;
    int32_t device_count = getDeviceCount();
    cudaDeviceProp device_properties;
    BOOST_LOG_TRIVIAL(info) << "Device(s) memeory check:" << '\n';
    for(int i=0; i<device_count; i++)
    {
        cudaSetDevice(i);
        cudaGetDeviceProperties(&device_properties, i);
        cudaMemGetInfo(&free, &total);
        BOOST_LOG_TRIVIAL(info)  << "  Device " << i+1 << ", " << device_properties.name  << ": " << (free>required_size_MB ? "OK" : "Not enough") << '\n';
        if(free<required_size_MB)
            BOOST_LOG_TRIVIAL(fatal) << "Not enough GPU memory for the simulation in device "<< i <<"! Required=" << required_size_MB <<" MB, Available=" << free << " MB";
        memory_ok = free<required_size_MB ? false:memory_ok;
    }
    return memory_ok;
}


#endif // _HELPER_H_
