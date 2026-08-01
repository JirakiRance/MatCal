#include <cmath>

double matcal_calculus_multi_tu_value();
double matcal_ode_multi_tu_value();

int main() {
    return std::abs(matcal_calculus_multi_tu_value() + matcal_ode_multi_tu_value() - 1.5) < 1.0e-12 ? 0 : 1;
}
