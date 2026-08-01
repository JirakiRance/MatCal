#include "Basics.hpp"
#include "Matrix.hpp"

double matcal_multi_tu_norm_b() {
    MatCal::Utils::Matrix m({{1.0, -2.0}, {3.0, 4.0}});
    return MatCal::Algorithm::Matrix::norm_infinite(m);
}
