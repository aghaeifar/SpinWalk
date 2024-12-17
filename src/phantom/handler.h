

#include <string>
#include <cstdint>

namespace phantom {
    struct execute_args {
        bool cylinder = false;
        bool sphere = false;
        float radius = 0.0f;
        float orientation = 0.0f;
        float volume_fraction = 0.0f;
        float fov = 0.0f;
        uint32_t resolution = 0;
        float dchi = 0.0f;
        float oxy_level = 0.0f;
        int32_t seed = 0;
        std::string output = "default.h5"; // Default value for name
    };

    class handler {
    public:
        static bool execute(const execute_args& args);
    };
}
