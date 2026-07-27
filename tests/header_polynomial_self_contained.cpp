#include "MatCal/Polynomial/Polynomial.hpp"

int matcal_polynomial_header_probe() {
    MatCal::Polynomial::Polynomial p{1.0, 2.0};
    return p.evaluate(2.0) == 5.0 ? 0 : 1;
}
