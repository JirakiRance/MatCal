#include "MatCal/LeastSquares/LeastSquares.hpp"

int header_least_squares_self_contained() {
    auto result = MatCal::LeastSquares::fit_polynomial_degree(0, {0.0, 1.0}, {2.0, 2.0});
    return result.success() ? 0 : 1;
}
