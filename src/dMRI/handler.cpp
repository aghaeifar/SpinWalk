
#include <iostream>

#include "handler.h"
#include "pgse.h"


namespace dMRI {
    bool handler::execute(const execute_args& args) {        
        std::cout << "Generating PGSE gradient table...\n";
        pgse gg;
        gg.set_parameters(args.b_value, args.start_ms, args.delta_ms, args.DELTA_ms);
        gg.set_direction(args.dir);
        if(gg.run(args.output)==false)
            return false;

        std::cout << "Diffusion gradient table is generated successfully." << '\n';
        return true;
    }
}