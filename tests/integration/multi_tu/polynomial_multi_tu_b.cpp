#include "MatCal/Polynomial/Polynomial.hpp"

double matcal_polynomial_multi_tu_b(const MatCal::Polynomial::Polynomial& p) {
    return p.derivative().evaluate(3.0);
}
