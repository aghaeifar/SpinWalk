cmake_minimum_required(VERSION 3.18)
message(STATUS "CMake version: ${CMAKE_VERSION}")
include(CheckLanguage)

set(project "SpinWalk")
project(${project})

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Find TBB
find_package(TBB REQUIRED)

# Boost
if(CMAKE_VERSION VERSION_GREATER "3.30")
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost COMPONENTS program_options log log_setup REQUIRED)

check_language(CUDA)
# Check for CUDA-capable GPU
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    # set_property(TARGET ${target} PROPERTY CUDA_ARCHITECTURES native)
    set(CMAKE_CUDA_ARCHITECTURES 75 86 89)
else()
    message("CUDA not found. Compiling without GPU support.")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add Boost and CUDA include directories
include_directories(${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ./include ./src ./src/shapes)
# Add the executable
set(SOURCES ./src/spinwalk.cu ./src/kernels.cu ./src/file_utils.cpp ./src/shapes/shape_base.cpp ./src/shapes/cylinder.cpp ./src/shapes/sphere.cpp)

if(CMAKE_CUDA_COMPILER)
    add_executable(${project} ${SOURCES})
    target_link_libraries(${project} ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_LIBRARIES} TBB::tbb)
else()    
    set(RENAMED_SOURCES)
    foreach(OLD_FILE ${SOURCES})
        # Get the file name without extension
        get_filename_component(DIR ${OLD_FILE} DIRECTORY)
        get_filename_component(FILE_NAME_WE ${OLD_FILE} NAME_WE)
        # Define the new file name
        set(NEW_FILE ${DIR}/${FILE_NAME_WE}.cpp)        
        # Copy and rename the file
        file(COPY_FILE ${OLD_FILE} ${NEW_FILE})
        list(APPEND RENAMED_SOURCES ${NEW_FILE})
    endforeach()
    add_executable(${project} ${RENAMED_SOURCES})
    target_link_libraries(${project} ${Boost_LIBRARIES} ${HDF5_CXX_LIBRARIES} TBB::tbb)
endif()

set_target_properties(${project} PROPERTIES OUTPUT_NAME "spinwalk")

if (UNIX)
    install(TARGETS ${project} DESTINATION bin) 
endif()

# cmake ..
# cmake --build . --config Release
