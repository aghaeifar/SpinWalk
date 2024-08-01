
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: file_utils.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

#ifndef __FILE_UTILS_H__
#define __FILE_UTILS_H__

#include <map>

#include "simulation_parameters.h"


//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
template <typename T>
size_t product(const std::vector<T>& v)
{
    size_t n_elements = v.size()==0 ? 0 : 1;
    for (const auto& e: v)
        n_elements *= e;
    return n_elements;
}


namespace file_utils
{
//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
bool read_config(std::string config_filename, simulation_parameters *param, std::vector<double>& sample_length_scales, std::map<std::string, std::vector<std::string> >& filenames, bool isParentConfig = false);

//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
bool read_phantom(std::string phantom_filename, std::vector<float> &fieldmap, std::vector<uint8_t> &mask, simulation_parameters *param);

//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
std::vector<size_t> get_size_h5(std::string input_filename, std::string dataset_name);

//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
bool read_h5(std::string input_filename, void *data, std::string dataset_name, std::string data_type);

//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
bool save_h5(std::string output_filename, void *data, std::vector<size_t> dims, std::string dataset_name, std::string data_type);

}
#endif  // __FILE_UTILS_H__
