#include "MatCal/Interpolation/CubicSpline.hpp"

int matcal_cubic_spline_header_self_contained() {
    MatCal::Interpolation::CubicSpline spline({0.0, 1.0}, {0.0, 2.0});
    return spline.evaluate(0.5) == 1.0 ? 0 : 1;
}
