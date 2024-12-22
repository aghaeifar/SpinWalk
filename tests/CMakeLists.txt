cmake_minimum_required(VERSION 3.24)
include(CTest)

set(project_test "SpinWalk_test")
project(${project_test})

# Boost
if(CMAKE_VERSION VERSION_GREATER "3.30")
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED COMPONENTS unit_test_framework log)

if (NOT Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message("Boost Unit Test Framework not found. Not compiling tests")
    return()
endif ()

if (NOT Boost_USE_STATIC_LIBS)
    add_definitions(-DBOOST_TEST_DYN_LINK)
endif()


# include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include ${Boost_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIRS})

add_executable(${project_test} test_phantom.cpp test_config.cpp)
target_link_libraries(${project_test} ${Boost_LIBRARIES} SpinWalk_lib)
set_target_properties(${project_test} PROPERTIES OUTPUT_NAME "spinwalk_test")

add_test(NAME check COMMAND spinwalk_test)

if (UNIX)
    install(TARGETS ${project_test} DESTINATION bin) 
endif()

# spinwalk_test --log_level=test_suite --report_level=short
