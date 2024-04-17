
#ifndef _HELPER_H_
#define _HELPER_H_


#include <iostream>
#include <boost/log/trivial.hpp>
#include "helper_cuda.h"


#define SPINWALK_VERSION_MAJOR 1
#define SPINWALK_VERSION_MINOR 9
#define SPINWALK_VERSION_PATCH 0

//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
inline void print_logo()
{ 
 std::cout << " \n"
" ____            _          __        __          _   _        \n"
"/ ___|   _ __   (_)  _ __   \\ \\      / /   __ _  | | | | __    \n"
"\\___ \\  | '_ \\  | | | '_ \\   \\ \\ /\\ / /   / _` | | | | |/ /    \n"
" ___) | | |_) | | | | | | |   \\ V  V /   | (_| | | | |   <     \n"
"|____/  | .__/  |_| |_| |_|    \\_/\\_/     \\__,_| |_| |_|\\_\\    \n"
"        |_|                                                    \n\n";

std::cout << "SpinWalk ver. " << SPINWALK_VERSION_MAJOR << "." << SPINWALK_VERSION_MINOR << "." << SPINWALK_VERSION_PATCH << std::endl;
}

//---------------------------------------------------------------------------------------------
//  check for CUDA and GPU device
//---------------------------------------------------------------------------------------------

uint32_t getDeviceCount()
{
    int32_t device_count;
    checkCudaErrors(cudaGetDeviceCount(&device_count));
    return device_count;
}

void print_device_info()
{
    const int kb = 1024;
    const int mb = kb * kb;
    size_t free, total;
    
    int32_t cuda_version, driver_version;
    int32_t device_count = getDeviceCount();
    cudaRuntimeGetVersion(&cuda_version);
    cudaDriverGetVersion(&driver_version);
    std::cout << "\nDriver version: "<< driver_version << ", CUDA version: "<< cuda_version << ", Number of devices: " << device_count << std::endl;

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
    size_t free, total;
    bool memory_ok = true;
    int32_t device_count = getDeviceCount();
    cudaDeviceProp device_properties;
    std::cout << "Device(s) memeory check:" << '\n';
    for(int i=0; i<device_count; i++)
    {
        cudaSetDevice(i);
        cudaGetDeviceProperties(&device_properties, i);
        cudaMemGetInfo(&free, &total);
        std::cout << "  Device " << i+1 << ", " << device_properties.name  << ": " << (free>required_size_MB ? "OK" : "Not enough") << '\n';
        if(free<required_size_MB)
            BOOST_LOG_TRIVIAL(fatal) << "Not enough GPU memory for the simulation in device "<< i <<"! Required=" << required_size_MB <<" MB, Available=" << free << " MB";
        memory_ok = free<required_size_MB ? false:memory_ok;
    }
    return memory_ok;
}


#endif // _HELPER_H_
