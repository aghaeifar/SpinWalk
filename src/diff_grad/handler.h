

#include <string>
#include <vector>

namespace diff_grad {
    struct execute_args {
        float b_value = 0.f;        
        uint32_t start_ms = 0;        
        uint32_t delta_ms = 0;
        uint32_t DELTA_ms = 0;
        std::vector<float> dir = {0.f, 0.f, 1.f};
        std::string output = "default.h5"; // Default value for name
    };

    class handler {
    public:
        static bool execute(const execute_args& args);
    };
}