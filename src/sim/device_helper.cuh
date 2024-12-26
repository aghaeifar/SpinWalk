
#ifndef DEVICE_HELPER_CUH
#define DEVICE_HELPER_CUH

#include <string>

namespace sim
{
bool check_CUDA();

uint32_t get_device_count();

std::string convert_cuda_version(int32_t value) ;

void print_device_info();

bool check_memory_size(size_t required_size_MB);

} // namespace sim
#endif // DEVICE_HELPER_CUH
