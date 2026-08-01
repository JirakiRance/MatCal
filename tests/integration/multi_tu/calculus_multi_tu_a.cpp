#include "MatCal/Calculus/Calculus.hpp"

double matcal_calculus_multi_tu_value() {
    return MatCal::Calculus::integrate_newton_cotes([](double x) { return x; }, 0.0, 1.0, 1).value;
}
