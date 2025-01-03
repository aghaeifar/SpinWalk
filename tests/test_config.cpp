#define BOOST_TEST_MODULE "SpinWlkTest"
// boost headers
#include <boost/test/unit_test.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

#include <filesystem>

#include "../src/config/config_generator.h"
#include "../include/ini.h"

struct InitLogging {
    InitLogging() { 
    std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "spinwalk_test";
    std::string log_filename = temp_dir / "spinwalk_test.log";
    auto fileSink = boost::log::add_file_log(boost::log::keywords::file_name=log_filename, boost::log::keywords::target_file_name = log_filename, boost::log::keywords::format = "[%TimeStamp%] [%Severity%]: %Message%", boost::log::keywords::auto_flush = true);
    boost::log::add_common_attributes();
    std::cout << "Log file location: " << log_filename << '\n';
    }
};

BOOST_GLOBAL_FIXTURE(InitLogging);

struct welcome_msg_config_module {
    welcome_msg_config_module() { 
        BOOST_LOG_TRIVIAL(info) << "--------- test_config_module ---------";
    }
};

BOOST_FIXTURE_TEST_SUITE(test_config_module, welcome_msg_config_module)

// BOOST_AUTO_TEST_SUITE(test_config_module)

BOOST_AUTO_TEST_CASE(config_creation) 
{
    boost::log::core::get()->set_logging_enabled(false);

    std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "spinwalk_test";
    std::filesystem::path temp_file = temp_dir / "gre.ini";
    std::filesystem::path temp_file_parent = temp_dir / "default_config.ini";
    std::vector<std::string> phantoms = {"phantom1.h5", "phantom2.h5"};
    uint32_t timestep_us = 25;
    uint32_t TE = 12345;

    config::config_generator cgen;
    BOOST_REQUIRE(cgen.generate_gre(TE, timestep_us, phantoms, temp_file.string()));
    BOOST_REQUIRE(std::filesystem::exists(temp_file));
    BOOST_REQUIRE(std::filesystem::exists(temp_file_parent));

    mINI::INIFile file(temp_file.string());
    mINI::INIStructure ini;
    BOOST_REQUIRE(file.read(ini));
    BOOST_CHECK_EQUAL(std::stoi(ini["SCAN_PARAMETERS"]["TE[0]"]), TE);
    BOOST_CHECK_EQUAL(std::stoi(ini["SCAN_PARAMETERS"]["TIME_STEP"]), timestep_us);

    std::filesystem::remove(temp_file);
    std::filesystem::remove(temp_file_parent);
    BOOST_CHECK(!std::filesystem::exists(temp_file));
    BOOST_CHECK(!std::filesystem::exists(temp_file_parent));
}

BOOST_AUTO_TEST_SUITE_END()