#include "Insert.hpp"
#include "MatCal/Interpolation/LinearInterpolator.hpp"
#include "MatCal/Polynomial/Polynomial.hpp"
#include "MatCal/Roots/Roots.hpp"

int main() {
    auto root = MatCal::Roots::solve_bisection([](double x) { return x - 2.0; }, 0.0, 4.0);
    if (!root.success()) {
        return 1;
    }

    MatCal::Interpolation::LinearInterpolator line({0.0, 1.0}, {1.0, 3.0});
    if (line.evaluate(0.5) != 2.0) {
        return 2;
    }

    MatCal::Polynomial::Polynomial p{1.0, 2.0};
    if (p.evaluate(3.0) != 7.0) {
        return 3;
    }

    MatCal::Algorithm::Insert::LinearInsert legacy({{0.0, 1.0}, {1.0, 3.0}});
    return legacy.calculate(0.5) == 2.0 ? 0 : 4;
}
