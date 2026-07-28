#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Interpolation/LinearInterpolator.hpp"
#include "MatCal/Polynomial/Polynomial.hpp"
#include "MatCal/Roots/Roots.hpp"
#include "QinJiuShao.hpp"

int main() {
    MatCal::Polynomial::Polynomial p{1.0, 2.0};
    if (p.evaluate(3.0) != 7.0) {
        return 1;
    }

    auto root = MatCal::Roots::solve_bisection([](double x) { return x - 2.0; }, 0.0, 4.0);
    if (!root.success() || root.value != 2.0) {
        return 5;
    }

    MatCal::Interpolation::LinearInterpolator line({0.0, 1.0}, {1.0, 3.0});
    if (line.evaluate(0.5) != 2.0) {
        return 6;
    }

    MatCal::Utils::QinJiuShao legacy({{1, 2.0}, {0, 1.0}});
    if (legacy.calculate(3.0) != 7.0) {
        return 2;
    }

    auto result = MatCal::Linalg::solve_dense_partial_pivot(
        MatCal::Linalg::DenseMatrix{{2.0, 0.0}, {0.0, 4.0}},
        MatCal::Linalg::Vector{6.0, 8.0});
    if (!result.success()) {
        return 3;
    }
    return result.solution[0] == 3.0 && result.solution[1] == 2.0 ? 0 : 4;
}
