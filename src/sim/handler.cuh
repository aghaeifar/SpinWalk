

#include <string>
#include <vector>
#include <cstdint>

namespace sim {
    struct execute_args {            
        bool use_cpu = false;        
        uint32_t device_id = 0;
        std::vector<std::string>  config_files; // Default value for name
    };

    class handler {
    public:
        static bool execute(const execute_args& args);
    };
}
