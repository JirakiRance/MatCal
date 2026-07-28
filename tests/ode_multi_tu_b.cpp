#include "MatCal/ODE/ODE.hpp"

double matcal_ode_multi_tu_value() {
    auto result = MatCal::ODE::rk4_step(
        [](double, const std::vector<double>& y) { return y; },
        0.0,
        std::vector<double>{1.0},
        0.0);
    return result.state[0];
}
