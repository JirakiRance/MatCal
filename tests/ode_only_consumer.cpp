#include "MatCal/ODE/ODE.hpp"

int main() {
    auto result = MatCal::ODE::rk4_step(
        [](double, const std::vector<double>& y) { return y; },
        0.0,
        std::vector<double>{1.0},
        0.1);
    return result.success() && result.state.size() == 1 ? 0 : 1;
}
