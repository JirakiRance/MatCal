#include "MatCal/LeastSquares/LeastSquares.hpp"

int least_squares_multi_tu_b() {
    auto result = MatCal::LeastSquares::fit_polynomial_degree(1, {0.0, 1.0}, {1.0, 3.0});
    return result.success() ? 0 : 1;
}
