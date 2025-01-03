/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: spinwalk.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : 
 * -------------------------------------------------------------------------- */
// standard libraries
#include <iomanip>
#include <filesystem>

// boost headers
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

// 3rd party headers
#include "CLI11.hpp"

// custom headers
#include "definitions.h"
#include "phantom/handler.h"
#include "dwi/handler.h"
#include "config/handler.h"
#include "sim/handler.cuh"
#include "sim/device_helper.cuh"

namespace bl = boost::log;

int main(int argc, char * argv[])
{
    // ========== command line arguments ==========
    bool arg_cyl = false, arg_sphere = false, use_cpu = false;
    std::string arg_output, config_file, arg_seqname;
    float arg_radius=50.f, arg_fov=1000.f, arg_dchi=0.11e-6, arg_oxy_level=0.75, arg_ori=90.f, arg_vol_fra=4.f;
    uint32_t arg_res=500, device_id=0, TE_us=1000, timestep_us=10;
    int32_t arg_seed = -1;
    std::vector<std::string> config_files, phantom_files;
    std::vector<uint32_t> bdelta;
    std::vector<float> bvector;
    std::vector<double> bvalue;
    
    // ========== parse command line arguments ==========
    CLI::App app{""};
    app.get_formatter()->column_width(40);
    app.get_formatter()->label("REQUIRED", "");
    app.set_version_flag("-v,--version", SPINWALK_VERSION);
    auto callback_gpu_info = [](int count){sim::print_device_info();  exit(0);};
    app.add_flag("-g,--gpu_info", callback_gpu_info, "Print GPU information");

    auto subcommand_sim = app.add_subcommand("sim", "Run Monte-Carlo simulation");
    subcommand_sim->add_option("-c,--configs", config_files, "Config. files as many as you want. e.g. -c config1.ini config2.ini ... configN.ini")->mandatory(true)->check(CLI::ExistingFile);
    subcommand_sim->add_flag("-p,--use_cpu", use_cpu, "Only run on CPU (default: GPU)");
    subcommand_sim->add_option("-d,--device", device_id, "Select GPU device (if there are multiple GPUs)");

    auto subcommand_phantom = app.add_subcommand("phantom", "Generate numerical phantom");
    subcommand_phantom->add_flag("-c,--cylinder", arg_cyl, "Fill phantom with cylinders");
    subcommand_phantom->add_flag("-s,--sphere", arg_sphere, "Fill phantom with spheres");
    subcommand_phantom->add_option("-r,--radius", arg_radius, "Radius of the cylinders/spheres in \u00B5m (negative value = random but smaller than radius)")->capture_default_str();
    subcommand_phantom->add_option("-n,--orientation", arg_ori, "Orientation of the cylinders in degree with respect to B0")->capture_default_str();
    subcommand_phantom->add_option("-v,--volume_fraction", arg_vol_fra, "Fraction of shapes volume to FoV volume in % <0.0 100.0>")->capture_default_str();
    subcommand_phantom->add_option("-f,--fov", arg_fov, "Voxel field of view in \u00B5m (isotropic)")->mandatory(true)->check(CLI::PositiveNumber);
    subcommand_phantom->add_option("-z,--resolution", arg_res, "Base resolution")->mandatory(true)->check(CLI::PositiveNumber);
    subcommand_phantom->add_option("-d,--dchi", arg_dchi, "Susceptibility difference between fully deoxygenated blood and tissue (default: 0.11e-6 in cgs units)")->capture_default_str();
    subcommand_phantom->add_option("-y,--oxy_level", arg_oxy_level, "Blood oxygenetation level <0.0 1.0> (-1 = exclude off-resonance effect and only generate the mask)")->capture_default_str();
    subcommand_phantom->add_option("-e,--seed", arg_seed, "Seed for random number generator in phantom creator (-1 = random seed)")->capture_default_str();
    subcommand_phantom->add_option("-o,--output", arg_output, "Path to save phantom (h5 format)")->mandatory(true);

    auto subcommand_config = app.add_subcommand("config", "Generate configuration file");
    subcommand_config->add_option("-s,--seq_name", arg_seqname, "Sequence name: GRE, SE, and bSSFP")->mandatory(true); 
    subcommand_config->add_option("-p,--phantoms", phantom_files, "Path to phantom files as many as you want. e.g. -p phantom1.h5 phantom2.h5 ... phantomN.h5")->mandatory(true); // must not check for existing file here, its path is relative to config file location
    subcommand_config->add_option("-e,--TE", TE_us, "Echo time in \u00B5s")->mandatory(true)->check(CLI::PositiveNumber);; 
    subcommand_config->add_option("-t,--timestep", timestep_us, "timestep in \u00B5s")->mandatory(true)->check(CLI::PositiveNumber);; 
    subcommand_config->add_option("-o,--output",config_file, "Path to save the configuration file")->mandatory(true);

    auto subcommand_diffusion = app.add_subcommand("dwi", "Generate diffusion gradient table");
    subcommand_diffusion->add_option("-b,--bvalue", bvalue, "b-value(s) as many as you want. e.g. -b 100 500 1000 5000 (s/mm\u00B2)")->mandatory(true);
    subcommand_diffusion->add_option("-v,--bvector", bvector, "Gradient direction: X Y Z, e.g. 0.267 0.534 0.801")->mandatory(true)->expected(3);
    subcommand_diffusion->add_option("-d,--delta", bdelta, "start time, \xCE\xB4 and \xCE\x94 in ms, e.g. 10 3 5 ")->mandatory(true)->expected(3);
    subcommand_diffusion->add_option("-c,--config", config_file, "input config file to insert PGSE gradients and excitation and refocusing RF")->mandatory(true)->check(CLI::ExistingFile);

    CLI11_PARSE(app, argc, argv);
    if(app.count_all() == 1){
        std::cout << app.help() << '\n';
        return 0;
    }
    if (subcommand_phantom->parsed() && arg_cyl == arg_sphere){  
        if (arg_cyl == true)      
            std::cout << "Error! select either --cylinder or --sphere, not both!"<< '\n';
        std::cout << subcommand_phantom->help() << '\n';
        return 0;
    }

    // ========== welcome ==========
    std::cout << SPINWALK_ASCII_ART << '\n';
    std::cout << "SpinWalk Version: " << SPINWALK_VERSION << '\n';  

    // ========== setup log ==========
    std::string log_filename = "spinwalk_" + std::to_string(device_id) + ".log";
    auto fileSink = bl::add_file_log(bl::keywords::file_name=log_filename, bl::keywords::target_file_name = log_filename, bl::keywords::format = "[%TimeStamp%] [%Severity%]: %Message%", bl::keywords::auto_flush = true);
    bl::add_common_attributes();
    std::cout << "Log file location: " << std::filesystem::current_path() / log_filename << '\n';

    // ========== generate diffusion gradients ==========
    if (subcommand_diffusion->parsed())
        if (dMRI::handler::execute({.start_ms=bdelta[0], .delta_ms=bdelta[1], .DELTA_ms=bdelta[2], .dir=bvector, .b_value=bvalue, .output=config_file}) == false){
            std::cout << ERR_MSG << "Diffusion gradient generation failed. See the log file " << log_filename <<", Aborting...!" << "\n";
            return 1;
        }

    // ========== generate config ==========
    if (subcommand_config->parsed())
        if (config::handler::execute({.seq_name=arg_seqname, .TE_us=TE_us, .timestep_us=timestep_us, .phantoms=phantom_files, .output=config_file}) == false){
            std::cout << ERR_MSG << "Configuration file generation failed. See the log file " << log_filename <<", Aborting...!" << "\n";
            return 1;
        }

    // ========== generate phantom ==========
    if(subcommand_phantom->parsed())
        if (phantom::handler::execute({.cylinder=arg_cyl, .sphere=arg_sphere, .radius=arg_radius, .orientation=arg_ori, .volume_fraction=arg_vol_fra, .fov=arg_fov, .resolution=arg_res, .dchi=arg_dchi, .oxy_level=arg_oxy_level, .seed=arg_seed, .output=arg_output}) == false){
            std::cout << ERR_MSG << "Phantom generation failed. See the log file " << log_filename <<", Aborting...!" << "\n";
            return 1;
        }
  
    // ========== Monte-Carlo simulation ==========
    if(subcommand_sim->parsed()){
         if (sim::handler::execute({.use_cpu=use_cpu, .device_id=device_id, .config_files=config_files}) == false){
            std::cout << ERR_MSG << "Simulation failed. See the log file " << log_filename <<", Aborting...!" << "\n";
            return 1;    
         }    
    }

    return 0;
}
