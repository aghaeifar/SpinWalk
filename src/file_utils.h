
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
#include <vector>
#include "miscellaneous.h"

namespace file_utils
{

//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------   
typedef struct output_header
{
    int32_t dim1, dim2, dim3, dim4;
    output_header(int32_t d1, int32_t d2=1, int32_t d3=1, int32_t d4=1): dim1(d1), dim2(d2), dim3(d3), dim4(d4){}
} output_header;


typedef struct input_header
{
    uint32_t fieldmap_size[3];
    double sample_length[3];
    bool fieldmap_exists=true, mask_exists=true;
    size_t file_size=0;
    input_header(uint32_t *a, double *b) {memcpy(fieldmap_size, a, 3*sizeof(uint32_t)); memcpy(sample_length, b, 3*sizeof(sample_length[0])); }
    input_header(){};
} input_header;


//---------------------------------------------------------------------------------------------
//  
//---------------------------------------------------------------------------------------------
bool read_config(std::string config_filename, simulation_parameters& param, std::vector<double>& sample_length_scales, std::map<std::string, std::vector<std::string> >& filenames);

bool read_header(std::string filename, input_header &hdr_in);

bool read_fieldmap(std::string fieldmap_filename, std::vector<float> &fieldmap, std::vector<uint8_t> &mask, simulation_parameters& param);

bool read_file(std::string filename, std::vector<float> &storage);

bool save_output(char *data, size_t bytes, std::string output_filename, output_header hdr, std::vector<double> &additional_hdr);

}
#endif  // __FILE_UTILS_H__