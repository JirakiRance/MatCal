#include "MatCal/Interpolation/LinearInterpolator.hpp"

int matcal_linear_interpolator_header_self_contained() {
    MatCal::Interpolation::LinearInterpolator line({0.0, 1.0}, {0.0, 1.0});
    return line.evaluate(0.5) == 0.5 ? 0 : 1;
}
