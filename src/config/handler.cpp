
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <cctype>
#include <filesystem>

#include "handler.h"
#include "config_generator.h"


namespace config {

    bool handler::execute(const execute_args& args) {   
        std::string seq_name = args.seq_name;
        std::string output = std::filesystem::weakly_canonical(std::filesystem::absolute(args.output)).string();
        std::transform(args.seq_name.begin(), args.seq_name.end(), seq_name.begin(), [](unsigned char c) { return std::tolower(c); });     

        config_generator cgen;
        std::unordered_map<std::string, std::function<int(uint32_t, uint32_t, std::vector<std::string>, std::string)>> function_map = {
            {"gre", std::bind(&config_generator::generate_gre, &cgen, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
            {"se", std::bind(&config_generator::generate_se, &cgen, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)},
            {"bssfp", std::bind(&config_generator::generate_bssfp, &cgen, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4)}
        };

        // Lookup and call the corresponding function
        auto it = function_map.find(seq_name);
        if (it != function_map.end()) {
            if(it->second(args.TE_us, args.timestep_us, args.phantoms, output) == false) // Call the function with the arguments
                return false;
        } else {
            std::cout << "Invalid function name!\n";
            return false;
        }

        std::cout << "Configuration file is generated in " << output <<'\n';
        return true;
    }
}
