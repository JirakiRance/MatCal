#include <cmath>

double matcal_roots_multi_tu_value();
double matcal_interpolation_multi_tu_value();

int main() {
    const double value = matcal_roots_multi_tu_value() + matcal_interpolation_multi_tu_value();
    return std::abs(value - 5.0) < 1.0e-12 ? 0 : 1;
}
