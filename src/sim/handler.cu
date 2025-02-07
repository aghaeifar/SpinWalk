
#include <iostream>
#include <filesystem>
#include "handler.cuh"
#include "monte_carlo.cuh"

namespace sim {
    bool handler::execute(const execute_args& args) {    
        sim::monte_carlo mc(args.use_cpu, args.device_id); 
        for(const auto& config_file : args.config_files){
            std::cout << "<" << std::filesystem::path(config_file).filename().string() << ">\n";       
            if(mc.run(config_file) == false)
                return false;                
        }  
        std::cout << "Simulation completed successfully. See the log file\n";
        return true;
    }
}