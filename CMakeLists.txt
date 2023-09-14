cmake_minimum_required(VERSION 3.18)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(spinwalk LANGUAGES CUDA CXX)

# FindCUDA module
find_package(CUDA REQUIRED)

# Check for CUDA-capable GPU
if(CUDA_FOUND)
  # Check for GPU device
  cuda_select_nvcc_arch_flags(ARCH_FLAGS Auto)
  if(NOT ARCH_FLAGS)
    message("No CUDA-capable GPU found")
  endif()
else()
  message("CUDA not found")
endif()

# Find OpenMP
find_package(OpenMP)

# Set CUDA architecture (change according to your GPU)
# set(CUDA_ARCH "-arch=sm_75" CACHE STRING "CUDA architecture")
set_property(TARGET ${target} PROPERTY CUDA_ARCHITECTURES native)


# Set compiler flags for OpenMP
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_STANDARD 17)


# Add CUDA include directories
include_directories(./include)

# Add CUDA and OpenMP flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CUDA_ARCH} -Xcompiler=${OpenMP_CXX_FLAGS}")

# Add the executable
add_executable(spinwalk ./src/spinwalk.cu ./src/kernels.cu)

# Link CUDA and OpenMP libraries
target_link_libraries(spinwalk ${CUDA_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

if (UNIX)
    install(TARGETS spinwalk DESTINATION bin)    
endif()  

# cmake ..
# cmake --build . --config Release
