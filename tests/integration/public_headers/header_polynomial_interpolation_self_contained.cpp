#include "MatCal/Interpolation/PolynomialInterpolation.hpp"

int header_polynomial_interpolation_self_contained() {
    auto polynomial = MatCal::Interpolation::interpolate_lagrange({{0.0, 1.0}, {1.0, 3.0}});
    return polynomial.evaluate(0.5) > 0.0 ? 0 : 1;
}
