
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: file_utils.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */


#include <vector>

const float B0[3] = {0.f, 0.f, 1.f};

class shape
{
    public:
        shape();
        shape(float fov_um, size_t resolution, float dChi, float Y, std::string filename);
        ~shape();
        void set_space(float fov_um, size_t resolution);
        void set_blood_parameters(float dChi, float Y); 
        void set_filename(std::string filename = "shape.h5");      
        virtual bool run(){return true;}; 
        virtual bool save();
        virtual bool create_grid();
        virtual void generate_shapes() = 0;
        virtual void generate_mask_fieldmap() = 0;
        virtual void print_info();

    protected:
        float *m_pGrid;
        float *m_pFieldmap;
        uint8_t *m_pMask;
        size_t m_resolution;
        float m_fov, m_dChi, m_Y;
        float m_BVF; // blood volume fraction
        std::string m_filename;

    private: 
};


class cylinder : public shape
{
    public:
    cylinder();
    cylinder(float fov_um, size_t resolution, float dChi, float Y, float radius_um = 50, size_t num_cylinders = 5, float orientation = -1.0f, std::string filename = "shape.h5");
    ~cylinder();

    virtual bool run();
    virtual void set_cylinder_parameters(float radius_um = 50, size_t num_cylinders = 5, float orientation = -1.0f);    
    virtual void generate_shapes();
    virtual void generate_mask_fieldmap();
    virtual void print_info();

    protected:

    private: 
        float *m_pCylinder_points;
        float *m_pCylinder_directions; // normalized vector
        float *m_pCylinder_radius;
        float m_radius, m_orientation;
        size_t m_num_cylinders;
};
