#include "MatCal/Interpolation/PolynomialInterpolation.hpp"

int main() {
    auto polynomial = MatCal::Interpolation::interpolate_hermite({0.0, 1.0}, {0.0, 1.0}, {0.0, 2.0});
    return polynomial.evaluate(0.5) > 0.0 ? 0 : 1;
}
