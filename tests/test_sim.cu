#include <boost/test/unit_test.hpp>
#include <boost/log/trivial.hpp> 
#include <filesystem>

#include "../src/config/config_generator.h"
#include "../src/phantom/phantom_cylinder.h"
#include "../src/sim/monte_carlo.cuh"


struct welcome_msg_sim_module {
    welcome_msg_sim_module() { 
        BOOST_LOG_TRIVIAL(info) << "--------- test_sim_module ---------";
    }
};

BOOST_FIXTURE_TEST_SUITE(test_sim_module, welcome_msg_sim_module)


// BOOST_AUTO_TEST_SUITE(test_sim_module)

BOOST_AUTO_TEST_CASE(sim_gre) {
    std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "spinwalk_test";
    std::filesystem::path temp_file_gre = temp_dir / "gre_sim.ini";
    std::filesystem::path temp_file_parent = temp_dir / "default_config.ini";
    std::vector<std::string> phantoms = {(temp_dir / "phantom1.h5").string(), (temp_dir / "phantom2.h5").string()};


    phantom::cylinder cyl1(600, 300, 0.11e-6, 0.85, -20, 5, 0, 45, phantoms[0]);
    phantom::cylinder cyl2(600, 300, 0.11e-6, 0.78, -20, 5, 0, 45, phantoms[1]);
    BOOST_TEST(cyl1.run(true));
    BOOST_TEST(cyl2.run(true));

    uint32_t timestep_us = 25;
    uint32_t TE = 20000;

    config::config_generator cgen;
    BOOST_TEST(cgen.generate_gre(TE, timestep_us, phantoms, ""));
    cgen.add_param("SIMULATION_PARAMETERS", "NUMBER_OF_SPINS", "1000");
    BOOST_TEST(cgen.write_ini( temp_file_gre.string()));

    sim::monte_carlo mc(true, 0); 
    BOOST_TEST(mc.run(temp_file_gre.string()));
}

BOOST_AUTO_TEST_SUITE_END()