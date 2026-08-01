#include "MatCal/Calculus/Calculus.hpp"

int main() {
    auto result = MatCal::Calculus::integrate_newton_cotes([](double x) { return x; }, 0.0, 1.0, 1);
    return result.success() && result.value > 0.49 && result.value < 0.51 ? 0 : 1;
}
