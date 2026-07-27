#ifndef MATCAL_LINALG_TEST_SUPPORT_HPP
#define MATCAL_LINALG_TEST_SUPPORT_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace matcal_linalg_test {

inline int failures = 0;

inline void fail(const std::string& message) {
    ++failures;
    std::cerr << "FAIL: " << message << '\n';
}

inline void expect_true(bool condition, const std::string& message) {
    if (!condition) {
        fail(message);
    }
}

inline void expect_near(double actual, double expected, double tolerance, const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        fail(message + " expected " + std::to_string(expected) + " got " + std::to_string(actual));
    }
}

inline void expect_relative(double actual, double expected, double relative_tolerance, const std::string& message) {
    double scale = std::max(std::abs(expected), 1.0);
    expect_near(actual, expected, relative_tolerance * scale, message);
}

template <typename F>
void expect_throw(F&& func, const std::string& message) {
    try {
        func();
    } catch (const std::exception&) {
        return;
    }
    fail(message + " expected exception");
}

inline int finish(const std::string& name) {
    if (failures != 0) {
        std::cerr << failures << " failure(s) in " << name << '\n';
        return 1;
    }
    std::cout << name << " passed\n";
    return 0;
}

} // namespace matcal_linalg_test

#endif
