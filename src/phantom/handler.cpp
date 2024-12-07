

#include "handler.h"
#include "phantom_sphere.h"
#include "phantom_cylinder.h"

namespace phantom {
    bool handler::execute(const execute_args& args) {        
        if (args.cylinder){
            cylinder cyl(args.fov, args.resolution, args.dchi, args.oxy_level, args.radius, args.volume_fraction, args.orientation, args.seed, args.output);
            return cyl.run();
        }
        if (args.sphere){
            sphere sph(args.fov, args.resolution, args.dchi, args.oxy_level, args.radius, args.volume_fraction, args.seed,args.output);
            return sph.run();
        }

        return true;
    }
}