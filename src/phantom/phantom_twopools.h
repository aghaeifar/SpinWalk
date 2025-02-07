
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_twopools.h
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 27.01.2025
 * Descrip  : 
 * -------------------------------------------------------------------------- */

#ifndef TWOPOOLS_H
#define TWOPOOLS_H

#include "phantom_base.h"

namespace phantom
{
class twopools : public phantom_base
{
    public:
    twopools();
    twopools(float fov_um, size_t resolution, std::string filename = "shape.h5");
    virtual ~twopools();

    virtual bool run(bool write_to_disk) override;
    virtual bool generate_mask_fieldmap();

    friend std::ostream& operator<<(std::ostream& os, const twopools& obj);

};

}

#endif // TWOPOOLS_H