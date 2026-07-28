#include "MatCal/ODE/ODE.hpp"

int matcal_ode_header_self_contained() {
    auto result = MatCal::ODE::rk4_step(
        [](double, const std::vector<double>& y) { return y; },
        0.0,
        std::vector<double>{1.0},
        0.1);
    return result.success() ? 0 : 1;
}
