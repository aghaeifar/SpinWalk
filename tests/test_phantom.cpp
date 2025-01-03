#include <boost/test/unit_test.hpp>
#include <boost/log/trivial.hpp> 
#include "../src/phantom/phantom_sphere.h"
#include "../src/phantom/phantom_cylinder.h"


struct welcome_msg_phantom_module {
    welcome_msg_phantom_module() { 
        BOOST_LOG_TRIVIAL(info) << "--------- test_phantom_module ---------";
    }
};

BOOST_FIXTURE_TEST_SUITE(test_phantom_module, welcome_msg_phantom_module)

BOOST_AUTO_TEST_CASE(cylinder_creation) {
    float diff_margin = 2.f;
    float vf = 10.0f;
    float angles = 10.0f;
    phantom::cylinder cyl(600, 300, 0.11e-6, -1, -20, vf, 0, angles, "");
    BOOST_TEST(cyl.run(false));
    BOOST_TEST(std::abs(cyl.get_actual_volume_fraction() - vf) < diff_margin); 
}

BOOST_AUTO_TEST_CASE(sphere_creation) {
    float diff_margin = 2.f;
    float vf = 12.0;
    phantom::sphere sph(600, 300, 0.11e-6, -1, -20, vf, 0, "");
    BOOST_CHECK(sph.run(false));
    BOOST_CHECK(std::abs(sph.get_actual_volume_fraction() - vf) < diff_margin); 
}


BOOST_AUTO_TEST_SUITE_END()