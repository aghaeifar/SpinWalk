cmake_minimum_required(VERSION 3.24)
message(STATUS "CMake version: ${CMAKE_VERSION}")
include(CheckLanguage)

set(project_exe "SpinWalk_exe")
set(project_lib "SpinWalk_lib")
project(${project_exe})

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Find TBB
find_package(TBB REQUIRED)

# Boost
if(CMAKE_VERSION VERSION_GREATER "3.30")
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED COMPONENTS log log_setup)

check_language(CUDA)
# Check for CUDA-capable GPU
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    # set_property(TARGET ${target} PROPERTY CUDA_ARCHITECTURES native)
    # set(CMAKE_CUDA_ARCHITECTURES all)
    set(CMAKE_CUDA_ARCHITECTURES 75 86 89)
else()
    message("CUDA not found. Compiling without GPU support.")
endif()

find_package(OpenMP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_definitions(MINI_CASE_SENSITIVE)

# Enable AddressSanitizer for debug builds
if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options(-fsanitize=address -g -O1)
    add_link_options(-fsanitize=address)
endif()

# Add Boost and CUDA include directories
include_directories(${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ./include ./src)
# Add modules directory
add_subdirectory(src/dwi)
add_subdirectory(src/phantom)
add_subdirectory(src/config)
add_subdirectory(src/sim)

 # Add the sources files
set(SOURCE_MODULES ${SUBCOMMAND_PHANTOM} ${SUBCOMMAND_SIM} ${SUBCOMMAND_DWI} ${SUBCOMMAND_CONFIG})
set(SOURCE_MAIN ./src/spinwalk.cu)

# change extension of the files if CUDA is not found
if(NOT  CMAKE_CUDA_COMPILER) 
    set(RENAMED_SOURCES)
    foreach(OLD_FILE ${SOURCE_MODULES})
        # Get the file name without extension
        get_filename_component(DIR ${OLD_FILE} DIRECTORY)
        get_filename_component(FILE_NAME_WE ${OLD_FILE} NAME_WE)
        # Define the new file name
        set(NEW_FILE ${DIR}/${FILE_NAME_WE}.cpp)        
        # Copy and rename the file
        file(COPY_FILE ${OLD_FILE} ${NEW_FILE})
        list(APPEND RENAMED_SOURCES ${NEW_FILE})
    endforeach()
    set(SOURCE_MODULES ${RENAMED_SOURCES})

    file(COPY_FILE ./src/spinwalk.cu ./src/spinwalk.cpp)
    set(SOURCE_MAIN ./src/spinwalk.cpp)
endif()

# Add the executable
add_library(${project_lib} ${SOURCE_MODULES})
add_executable(${project_exe} ${SOURCE_MAIN})

# Add the libraries
if(CMAKE_CUDA_COMPILER)
    target_link_libraries(${project_lib} ${CUDA_LIBRARIES})
endif()
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP, adding to target link libraries.")
    target_link_libraries(${project_lib} OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP not found, skipping.")
endif()
target_link_libraries(${project_lib} ${Boost_LIBRARIES} ${HDF5_CXX_LIBRARIES} TBB::tbb)
target_link_libraries(${project_exe} PRIVATE ${project_lib})

# Set the output name
set_target_properties(${project_exe} PROPERTIES OUTPUT_NAME "spinwalk")
set_target_properties(${project_lib} PROPERTIES OUTPUT_NAME "spinwalk")


# Install the executable
if (UNIX)
    install(TARGETS ${project_exe} DESTINATION bin) 
    install(TARGETS ${project_lib} DESTINATION lib) 
endif()

# Add the tests directory
add_subdirectory(tests)

# cmake ..
# cmake --build . --config Release
