// --------------------------------------------------------------------------------------------
// generate default configuration file
// --------------------------------------------------------------------------------------------


#include <string>
#include <vector>

namespace config
{
bool generate_default_config(uint32_t TE_us = 10000, uint32_t timestep_us = 50, std::vector<std::string> phantoms = std::vector<std::string>(), std::string output="default_config.ini");
bool generate_gre(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output);
bool generate_se(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output);
bool generate_bssfp(uint32_t TE_us, uint32_t timestep_us, std::vector<std::string> phantoms, std::string output);

} // namespace config