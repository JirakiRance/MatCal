#include "MatCal/Polynomial/Polynomial.hpp"

MatCal::Polynomial::Polynomial matcal_polynomial_multi_tu_a();
double matcal_polynomial_multi_tu_b(const MatCal::Polynomial::Polynomial&);

int main() {
    auto p = matcal_polynomial_multi_tu_a();
    if (p.evaluate(3.0) != 7.0) {
        return 1;
    }
    return matcal_polynomial_multi_tu_b(p) == 6.0 ? 0 : 1;
}
