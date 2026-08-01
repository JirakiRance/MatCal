#include "MatCal/Polynomial/Polynomial.hpp"

MatCal::Polynomial::Polynomial matcal_polynomial_multi_tu_a() {
    return MatCal::Polynomial::Polynomial::from_terms({{2, 1.0}, {0, -2.0}});
}
