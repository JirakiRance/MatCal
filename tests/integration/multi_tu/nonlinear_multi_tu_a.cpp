#include "MatCal/Nonlinear/Nonlinear.hpp"

#include <vector>

int nonlinear_multi_tu_a() {
    auto result = MatCal::Nonlinear::solve_newton_system(
        [](const std::vector<double>& x) { return std::vector<double>{x[0] - 1.0}; },
        {0.0},
        [](const std::vector<double>&) {
            MatCal::Linalg::DenseMatrix j(1, 1);
            j(0, 0) = 1.0;
            return j;
        });
    return result.success() ? 0 : 1;
}
