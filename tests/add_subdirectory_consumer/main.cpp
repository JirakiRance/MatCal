#include <vector>

#include "Insert.hpp"
#include "MatCal/Calculus/Calculus.hpp"
#include "MatCal/Linalg/EigenSolvers.hpp"
#include "MatCal/Linalg/IterativeSolvers.hpp"
#include "MatCal/Interpolation/LinearInterpolator.hpp"
#include "MatCal/Interpolation/PolynomialInterpolation.hpp"
#include "MatCal/LeastSquares/LeastSquares.hpp"
#include "MatCal/Nonlinear/Nonlinear.hpp"
#include "MatCal/ODE/ODE.hpp"
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

    auto lagrange = MatCal::Interpolation::interpolate_lagrange({{0.0, 1.0}, {1.0, 3.0}});
    if (lagrange.evaluate(0.5) != 2.0) {
        return 7;
    }

    auto nonlinear = MatCal::Nonlinear::solve_newton_system(
        [](const std::vector<double>& x) { return std::vector<double>{x[0] - 2.0}; },
        {0.0},
        [](const std::vector<double>&) {
            MatCal::Linalg::DenseMatrix j(1, 1);
            j(0, 0) = 1.0;
            return j;
        });
    if (!nonlinear.success()) {
        return 8;
    }

    auto fit = MatCal::LeastSquares::fit_polynomial_degree(1, {0.0, 1.0}, {1.0, 3.0});
    if (!fit.success()) {
        return 9;
    }

    auto iterative = MatCal::Linalg::solve_gauss_seidel(
        MatCal::Linalg::DenseMatrix{{4.0, 1.0}, {2.0, 3.0}},
        MatCal::Linalg::Vector{1.0, 2.0});
    if (!iterative.success()) {
        return 10;
    }

    auto eigen = MatCal::Linalg::dominant_eigenpair(MatCal::Linalg::DenseMatrix::identity(1));
    if (!eigen.success()) {
        return 11;
    }

    MatCal::Polynomial::Polynomial p{1.0, 2.0};
    if (p.evaluate(3.0) != 7.0) {
        return 3;
    }

    auto integral = MatCal::Calculus::integrate_newton_cotes([](double x) { return x; }, 0.0, 1.0, 1);
    if (!integral.success()) {
        return 5;
    }

    auto step = MatCal::ODE::rk4_step(
        [](double, const std::vector<double>& y) { return y; },
        0.0,
        std::vector<double>{1.0},
        0.0);
    if (!step.success()) {
        return 6;
    }

    MatCal::Algorithm::Insert::LinearInsert legacy({{0.0, 1.0}, {1.0, 3.0}});
    return legacy.calculate(0.5) == 2.0 ? 0 : 4;
}
