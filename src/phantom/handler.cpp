

#include "handler.h"
#include "phantom_sphere.h"
#include "phantom_cylinder.h"

namespace phantom {
    bool handler::execute(const execute_args& args) {    
        bool status = true;    
        if (args.cylinder){
            std::cout << "Generating cylinder phantom..." << std::endl;
            cylinder cyl(args.fov, args.resolution, args.dchi, args.oxy_level, args.radius, args.volume_fraction, args.orientation, args.seed, args.output);
            status = status && cyl.run();
        }
        if (args.sphere){
            std::cout << "Generating sphere phantom..." << std::endl;
            sphere sph(args.fov, args.resolution, args.dchi, args.oxy_level, args.radius, args.volume_fraction, args.seed,args.output);
            status = status && sph.run();
        }
        std::cout << "Done." << std::endl;
        return true;
    }
}