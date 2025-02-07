

#include "handler.h"
#include "phantom_ply.h"
#include "phantom_sphere.h"
#include "phantom_cylinder.h"
#include "phantom_twopools.h"

namespace phantom {
    bool handler::execute(const execute_args& args) {    
        bool status = true;    
        if (args.cylinder){
            std::cout << "Generating cylinder phantom..." << std::endl;
            cylinder cyl(args.fov, args.resolution, args.dchi, args.oxy_level, args.radius, args.volume_fraction, args.orientation, args.seed, args.output);
            status = status && cyl.run(true);
        }
        if (args.sphere){
            std::cout << "Generating sphere phantom..." << std::endl;
            sphere sph(args.fov, args.resolution, args.dchi, args.oxy_level, args.radius, args.volume_fraction, args.seed,args.output);
            status = status && sph.run(true);
        }
        if(args.twopools){
            std::cout << "Generating phantom with two pools..." << std::endl;
            twopools tp(args.fov, args.resolution, args.output);
            status = status && tp.run(true);
        }
        if(args.ply){
            std::cout << "Generating phantom from triangular mesh..." << std::endl;
            ply ply_p(args.fov, args.resolution, args.dchi, args.oxy_level, args.ply_file, args.output);
            status = status && ply_p.run(true);
        }
        
        std::cout << "Done." << std::endl;
        return status;
    }
}