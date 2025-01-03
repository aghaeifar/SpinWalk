#include <cmath>

#include <boost/test/unit_test.hpp>
#include <boost/log/trivial.hpp> 

#include "../src/sim/kernels.cuh"


struct welcome_msg_test_kernel {
    welcome_msg_test_kernel() { 
        BOOST_LOG_TRIVIAL(info) << "--------- test_kernel ---------";
    }
};

BOOST_FIXTURE_TEST_SUITE(test_kernel, welcome_msg_test_kernel)

BOOST_AUTO_TEST_CASE(test_sub2ind_3d_row_major) {
    int64_t x = 1, y = 2, z = 3;
    int64_t len_dim_x = 10, len_dim_y = 10, len_dim_z = 10;

    int64_t expected = x * len_dim_z * len_dim_y + y * len_dim_z + z;
    BOOST_CHECK_EQUAL(sim::sub2ind(x, y, z, len_dim_x, len_dim_y, len_dim_z), expected);
}

BOOST_AUTO_TEST_CASE(test_sub2ind_4d_row_major) {
    int64_t x = 1, y = 2, z = 3, o = 4;
    int64_t len_dim_x = 10, len_dim_y = 10, len_dim_z = 10, len_dim_o = 10;

    int64_t expected = x * len_dim_o * len_dim_z * len_dim_y +
                       y * len_dim_z * len_dim_o +
                       z * len_dim_o +
                       o;
    BOOST_CHECK_EQUAL(sim::sub2ind(x, y, z, o, len_dim_x, len_dim_y, len_dim_z, len_dim_o), expected);
}

BOOST_AUTO_TEST_CASE(test_xrot) {
    float m0[3] = {0.0f, 0.0f, 1.0f};
    float m1[3] = {0.0f, 0.0f, 0.0f};
    float theta = 90.0f;
    float sin_theta = std::sin(theta * M_PI / 180.0f);
    float cos_theta = std::cos(theta * M_PI / 180.0f);

    sim::xrot(sin_theta, cos_theta, m0, m1);

    BOOST_TEST(std::abs(m1[0] - 0.0f) < 1e-5);
    BOOST_TEST(std::abs(m1[1] - (-1.0f)) < 1e-5);
    BOOST_TEST(std::abs(m1[2] - 0.0f) < 1e-5);
}

BOOST_AUTO_TEST_CASE(test_yrot) {
    float m0[3] = {1.0f, 0.0f, 0.0f};
    float m1[3] = {0.0f, 0.0f, 0.0f};
    float theta = 90.0f;
    float sin_theta = std::sin(theta * M_PI / 180.0f);
    float cos_theta = std::cos(theta * M_PI / 180.0f);

    sim::yrot(sin_theta, cos_theta, m0, m1);

    BOOST_TEST(std::abs(m1[0] - 0.0f) < 1e-5);
    BOOST_TEST(std::abs(m1[1] - 0.0f) < 1e-5);
    BOOST_TEST(std::abs(m1[2] - (-1.0f)) < 1e-5);
}

BOOST_AUTO_TEST_CASE(test_zrot) {
    float m0[3] = {1.0f, 0.0f, 0.0f};
    float m1[3] = {0.0f, 0.0f, 0.0f};
    float theta = 90.0f;
    float sin_theta = std::sin(theta * M_PI / 180.0f);
    float cos_theta = std::cos(theta * M_PI / 180.0f);

    sim::zrot(sin_theta, cos_theta, m0, m1);

    BOOST_TEST(std::abs(m1[0] - 0.0f) <  1e-5);
    BOOST_TEST(std::abs(m1[1] - 1.0f) <  1e-5);
    BOOST_TEST(std::abs(m1[2] - 0.0f) <  1e-5);
}

BOOST_AUTO_TEST_CASE(test_relax) {
    float m0[3] = {1.0f, 0.5f, -0.5f};
    float m1[3] = {0.0f, 0.0f, 0.0f};
    float e1 = 0.9f;
    float e2 = 0.8f;

    sim::relax(e1, e2, m0, m1);

    BOOST_TEST(std::abs(m1[0] - m0[0] * e2) <  1e-5);
    BOOST_TEST(std::abs(m1[1] - m0[1] * e2) <  1e-5);
    BOOST_TEST(std::abs(m1[2] - (1.0f + e1 * (m0[2] - 1.0f))) <  1e-5);
}

BOOST_AUTO_TEST_SUITE_END()