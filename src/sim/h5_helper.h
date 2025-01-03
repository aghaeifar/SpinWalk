
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: h5_helper.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : 
 * -------------------------------------------------------------------------- */

#ifndef __H5_HELPER_H__
#define __H5_HELPER_H__

#include <sstream>
#include <vector>
#include <fstream>
#include <highfive/highfive.hpp>
#include <boost/log/trivial.hpp> 

#include "simulation_parameters.cuh"

namespace sim
{
    template <typename T>
    size_t product(const std::vector<T>& v)
    {
        size_t n_elements = v.size()==0 ? 0 : 1;
        for (const auto& e: v)
            n_elements *= e;
        return n_elements;
    }

    class h5_helper
    {
    public:
        h5_helper() {};
        virtual ~h5_helper() {};
        template <typename T>
        static bool read(std::string filename, std::string dataset_name, bool resize_allowed, std::vector<T> &data);
        template <typename T>
        static bool write(std::string filename, std::string dataset_name, std::vector<size_t> dims, const std::vector<T> &data);
        static bool size(std::string filename, std::string dataset_name, std::vector<size_t> &dims);
    protected:
        static bool has_write_access(const std::string& path);
    };


    template <typename T>
    bool h5_helper::read(std::string filename, std::string dataset_name, bool resize_allowed, std::vector<T> &data)
    {
        BOOST_LOG_TRIVIAL(info) << "Reading dataset \"" << dataset_name << "\" from file \"" << filename << "\"";
        HighFive::File h5file(filename, HighFive::File::ReadOnly);
        if (h5file.exist(dataset_name) == false)
        {
            BOOST_LOG_TRIVIAL(error) << "dataset \"" << dataset_name << "\" does not exist in " << filename;
            return false;
        }

        HighFive::DataSet dataset = h5file.getDataSet(dataset_name);
        size_t new_size = product(dataset.getDimensions());
        if (new_size != data.size())
        {
            if (resize_allowed == false){            
                BOOST_LOG_TRIVIAL(error) << "dataset \"" << dataset_name << "\" has different size " << new_size << " vs " << data.size();
                return false;
            }
            data.resize(new_size);
        }        

        dataset.read_raw<T>(data.data());
        return true;
    }


    template <typename T>
    bool h5_helper::write(std::string filename, std::string dataset_name, std::vector<size_t> dims, const std::vector<T> &data)
    {
        std::ostringstream oss;
        for (const auto& elem : dims)
            oss << elem << " ";

        BOOST_LOG_TRIVIAL(info) << "Saving " << dataset_name << " with size = [" << oss.str() << "] and " << data.size() << " elements to: " << std::filesystem::absolute(filename);
        if(data.size() != product(dims))
        {
            BOOST_LOG_TRIVIAL(error) << "data size does not match the size of the dataset: " << data.size() << " vs " << product(dims);
            return false;
        }

        std::filesystem::path parent_path = std::filesystem::absolute(filename).parent_path();
        if (std::filesystem::is_directory(parent_path) == false)
        {
            BOOST_LOG_TRIVIAL(warning) << "cannot find directory " << parent_path.string() << ". Trying to create it.";
            if(std::filesystem::create_directories(parent_path) == false)
            {
                BOOST_LOG_TRIVIAL(error) << "cannot create directory " << parent_path.string();
                return false;
            }
        }

        if(has_write_access(filename) == false) {
            BOOST_LOG_TRIVIAL(error) << "write permission error for the file " << filename;
            return false;
        }

        // if file exists, open it in read-write mode, otherwise (re-)create it
        auto file_mode = std::filesystem::exists(filename) ? HighFive::File::ReadWrite : HighFive::File::Truncate;
        HighFive::File file(filename, file_mode);
        // if dataset exists, delete it
        if (file.exist(dataset_name))
            file.unlink(dataset_name);

        HighFive::DataSet dataset = file.createDataSet<T>(dataset_name, HighFive::DataSpace(dims));
        dataset.write_raw(data.data());
        return true;
    }


    bool h5_helper::size(std::string filename, std::string dataset_name, std::vector<size_t> &dims)
    {
        HighFive::File file(filename, HighFive::File::ReadOnly);
        if (file.exist(dataset_name) == false)
        {
            BOOST_LOG_TRIVIAL(error) << "dataset \"" << dataset_name << "\" does not exist in " << filename;
            return false;
        }
        HighFive::DataSet dataset = file.getDataSet(dataset_name);
        dims = dataset.getDimensions();
        return true;    
    }

    bool h5_helper::has_write_access(const std::string& path) {
        try {
            std::ofstream out(path, std::ios::app);
            if (out) {
                out.close();
                return true;
            }            
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << '\n';
        }
        return false;
    }

}

#endif  // __H5_HELPER_H__
