

#include <string>
#include <vector>

namespace config {
    struct execute_args {
        std::string seq_name = "default";        
        uint32_t TE_us = 0;
        uint32_t timestep_us = 0;
        std::vector<std::string> phantoms;   
        std::string output = "config_default.ini"; // Default value for name
    };

    class handler {
    public:
        static bool execute(const execute_args& args);
    };
}