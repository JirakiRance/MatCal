#include "MatCal/Roots/Roots.hpp"

double matcal_roots_multi_tu_value() {
    auto result = MatCal::Roots::solve_bisection([](double x) { return x - 2.0; }, 0.0, 4.0);
    return result.value;
}
