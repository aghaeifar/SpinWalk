
/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: miscellaneous.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#ifndef HOST_DEVICE_ARRAY_CUH
#define HOST_DEVICE_ARRAY_CUH

#include <vector>

#ifdef __CUDACC__
#include "helper_cuda.h"
#endif

#define CPU false
#define GPU true

template <typename T>
struct host_device_array {
    T *ptr, *d_ptr;
    size_t size;
    std::vector<T> data;

    host_device_array() : ptr(nullptr), d_ptr(nullptr), size(0){}
    ~host_device_array(){
#ifdef __CUDACC__
        if (d_ptr != nullptr)
            checkCudaErrors(cudaFree(d_ptr));
 #endif
        ptr = nullptr;
        d_ptr = nullptr;
    }

    void init(bool device = CPU, bool copy_from_host = false){
        bool resize = false;
        if (resize = size != data.size())
            size = data.size();
        if(device == CPU)
            ptr = data.data();
#ifdef __CUDACC__      
        if(device == GPU){
            if (resize && d_ptr != nullptr){
                checkCudaErrors(cudaFree(d_ptr));
                d_ptr = nullptr;
            }
            if (d_ptr == nullptr)
                checkCudaErrors(cudaMalloc(&d_ptr, data.size()*sizeof(T)));
            if(copy_from_host)
                checkCudaErrors(cudaMemcpy(d_ptr, data.data(), data.size()*sizeof(T), cudaMemcpyHostToDevice));
            ptr = d_ptr;
        }
#endif
    }

    void copy_to_host(){
#ifdef __CUDACC__
        if (d_ptr == nullptr)
            throw std::runtime_error("device pointer is not initialized");
        checkCudaErrors(cudaMemcpy(data.data(), d_ptr, data.size()*sizeof(T), cudaMemcpyDeviceToHost));
#endif
    }

};

#endif // HOST_DEVICE_ARRAY_CUH