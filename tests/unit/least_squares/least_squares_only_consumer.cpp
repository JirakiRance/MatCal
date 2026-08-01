#include "MatCal/LeastSquares/LeastSquares.hpp"

int main() {
    auto result = MatCal::LeastSquares::fit_polynomial_degree(1, {0.0, 1.0}, {1.0, 2.0});
    return result.success() ? 0 : 1;
}
