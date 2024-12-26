set(SUBCOMMAND_SIM
    ${CMAKE_CURRENT_SOURCE_DIR}/monte_carlo.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/handler.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/device_helper.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/config_reader.cpp
    PARENT_SCOPE)