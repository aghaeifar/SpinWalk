
#include <iostream>
#include <algorithm>

#include "handler.h"
#include "config_generator.h"


namespace config {

    bool handler::execute(const execute_args& args) {   
        std::string seq_name;
        std::transform(args.seq_name.begin(), args.seq_name.end(), seq_name.begin(), [](unsigned char c) { return std::tolower(c); });     

        if (seq_name == "default"){
            if(generate_default_config(args.output) == false)
                return false;
        } else {
            std::cout << "Not implemented yet\n";
            return false;
        }

        std::cout << "Default configuration file is generated in " << args.output <<'\n';
        return true;
    }
}
