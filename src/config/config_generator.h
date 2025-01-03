
#include <string>
#include <vector>
#include "ini.h"

namespace config
{
class config_generator{
    public:
        config_generator(){};
        ~config_generator(){};
        bool generate_default_config(uint32_t TE_us = 10000, uint32_t timestep_us = 50, std::vector<std::string> phantoms = std::vector<std::string>());
        bool generate_gre(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output);
        bool generate_se(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output);
        bool generate_bssfp(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output);
        void add_param(std::string section, std::string key, std::string value);
        bool write_ini(std::string output);
        
    private:
        mINI::INIStructure ini, ini_parent;
};

} // namespace config