#include <cmath>
#include <iostream>

double matcal_multi_tu_norm_a();
double matcal_multi_tu_norm_b();

int main() {
    const double a = matcal_multi_tu_norm_a();
    const double b = matcal_multi_tu_norm_b();
    if (std::abs(a - 6.0) > 1e-12 || std::abs(b - 7.0) > 1e-12) {
        std::cerr << "Unexpected multi-TU norm results: " << a << ", " << b << '\n';
        return 1;
    }
    std::cout << "MatCal multi-translation-unit link test passed\n";
    return 0;
}
