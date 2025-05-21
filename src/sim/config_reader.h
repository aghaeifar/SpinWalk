
#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include <string>
#include <vector>
#include <map>

struct parameters;
struct parameters_hvec;

namespace sim {
    enum e_scale_type {s_fov=0, s_gradient=1, s_phase_cycling=2};
    class config_reader {
    public:
        config_reader() {};
        virtual ~config_reader() {};
        virtual bool prepare(std::string config_filename, parameters *param, parameters_hvec *param_hvec);
        std::vector<std::string> get_filename(std::string key) const { return files_container.at(key); }        
        std::vector<float> get_scales() const { return scales; }
        e_scale_type get_scale_type() const { return scale_type; }
        std::string get_output_filename(size_t ind) const { return output_files.at(ind); }
        void dump() const;

    protected:
        bool read(std::string config_filename);
        bool check();
        void timing_scale();
        void cleanup();

    private:
        std::string config_filename = "";
        std::string output_dir = "";
        std::string seq_name = "";
        std::map<std::string, std::vector<std::string> > files_container = {{"PHANTOM", std::vector<std::string>()}, // input:  map of off-resonance in Tesla
                                                                            {"XYZ0", 	std::vector<std::string>()}, // input:  spins starting spatial positions in meters
                                                                            {"M0", 		std::vector<std::string>()}  // input:  spins initial magnetization
                                                                            };
        std::vector<std::string> output_files; // output: spins final magnetization + spatial positions in meters + tissue index
        std::vector<float> scales; 
        e_scale_type scale_type;
        parameters *param = nullptr;
        parameters_hvec *param_hvec = nullptr;
    };

} // namespace sim

#endif // CONFIG_READER_H