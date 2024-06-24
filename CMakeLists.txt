cmake_minimum_required(VERSION 3.18)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(spinwalk LANGUAGES CUDA CXX)

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)
# Find OpenMP
find_package(OpenMP REQUIRED)
# Boost
find_package(Boost COMPONENTS program_options log log_setup REQUIRED)
# Find TBB
find_package(TBB REQUIRED)
# FindCUDA module
find_package(CUDA)

# Check for CUDA-capable GPU
if(CUDA_FOUND)
    # Check for GPU device
    cuda_select_nvcc_arch_flags(ARCH_FLAGS Auto)
    if(NOT ARCH_FLAGS)
        message("No CUDA-capable GPU found")
    endif()
    # Set CUDA architecture (change according to your GPU)
    # set(CUDA_ARCH "-arch=sm_75" CACHE STRING "CUDA architecture")
    # set_property(TARGET ${target} PROPERTY CUDA_ARCHITECTURES native)
    set_property(TARGET ${target} PROPERTY CUDA_ARCHITECTURES 70 75 80 86 90)
    # Add CUDA and OpenMP flags
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CUDA_ARCH} -Xcompiler=${OpenMP_CXX_FLAGS}")
else()
    message("CUDA not found")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-lboost_program_options")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(Boost_USE_STATIC_LIBS OFF)

# Add Boost and CUDA include directories
include_directories(${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ./include ./src ./src/shapes)
# Add the executable
add_executable(spinwalk ./src/spinwalk.cu ./src/kernels.cu ./src/file_utils.cpp ./src/shapes/shape_base.cpp ./src/shapes/cylinder.cpp ./src/shapes/sphere.cpp)
# Link CUDA and OpenMP libraries
target_link_libraries(spinwalk ${CUDA_LIBRARIES} ${OpenMP_CXX_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_LIBRARIES} TBB::tbb)

if (UNIX)
    install(TARGETS spinwalk DESTINATION bin)    
endif()  

# cmake ..
# cmake --build . --config Release
