#ifndef MATCAL_TEST_SUPPORT_HPP
#define MATCAL_TEST_SUPPORT_HPP

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "Matrix.hpp"

namespace matcal_test {

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

template <typename F>
void expect_throw(F&& func, const std::string& message) {
    try {
        func();
    } catch (const std::exception&) {
        return;
    }
    fail(message + " expected exception");
}

inline MatCal::Utils::Matrix as_matrix(std::unique_ptr<MatCal::Utils::AbstractMatrix> value) {
    auto* matrix = dynamic_cast<MatCal::Utils::Matrix*>(value.get());
    if (!matrix) {
        throw std::runtime_error("Expected legacy operation to return MatCal::Utils::Matrix");
    }
    return *matrix;
}

inline int finish(const std::string& name) {
    if (failures != 0) {
        std::cerr << failures << " failure(s) in " << name << '\n';
        return 1;
    }
    std::cout << name << " passed\n";
    return 0;
}

} // namespace matcal_test

#endif
