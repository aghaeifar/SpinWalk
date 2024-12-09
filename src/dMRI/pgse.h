

#include <string>
#include <vector>

namespace dMRI
{

class pgse
{
public:
    pgse(){};
    virtual ~pgse(){};
    virtual void set_parameters(double b_value, uint32_t start_ms, uint32_t delta_ms, uint32_t DELTA_ms);
    virtual void set_direction(std::vector<float> dir){this->dir = dir;}
    virtual bool run(std::string config_file);

protected:
    bool read_timestep(std::string config_file, int &timestep_us);
    int read_RF(std::string config_file){return 0;} 

private:
    std::vector<float> dir;
    uint32_t DELTA_ms = 0, delta_ms = 0, start_ms = 0;
    double b_value = 0;
};

}