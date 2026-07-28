#include "MatCal/Interpolation/CubicSpline.hpp"
#include "MatCal/Interpolation/LinearInterpolator.hpp"

double matcal_interpolation_multi_tu_value() {
    MatCal::Interpolation::LinearInterpolator line({0.0, 1.0}, {1.0, 3.0});
    MatCal::Interpolation::CubicSpline spline({0.0, 1.0}, {0.0, 2.0});
    return line.evaluate(0.5) + spline.evaluate(0.5);
}
