#include "MatCal/Interpolation/PolynomialInterpolation.hpp"

int interpolation_remaining_multi_tu_c() {
    auto result = MatCal::Interpolation::interpolate_newton_finite(1.0, 0.0, {0.0, 1.0});
    return result.polynomial.evaluate(0.5) > 0.0 ? 0 : 1;
}
