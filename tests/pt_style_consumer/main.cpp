#include <vector>

#include "Insert.hpp"
#include "Iteration.hpp"

int main() {
    std::vector<double> xs{0.0, 1.0, 2.0};
    std::vector<double> ys{1.0, 3.0, 5.0};

    MatCal::Algorithm::Insert::LinearInsert linear(xs, ys);
    if (linear.calculate(0.5) != 2.0) {
        return 1;
    }

    MatCal::Algorithm::Insert::CubicSpline spline(xs, ys);
    if (spline.getXs().size() != 3 || spline.getYs().size() != 3 || spline.getM().size() != 3) {
        return 2;
    }
    if (spline.calculate(1.5) != 4.0) {
        return 3;
    }

    std::vector<std::pair<double, double>> interpolation_data{{0.0, 1.0}, {1.0, 2.0}, {2.0, 5.0}};
    MatCal::Algorithm::Insert::LagrangeInsert lagrange(interpolation_data);
    MatCal::Algorithm::Insert::NewtonInsert_Quotient quotient(interpolation_data);
    std::vector<double> finite_y{0.0, 1.0, 4.0};
    MatCal::Algorithm::Insert::NewtonInsert_Finite finite(1.0, 0.0, finite_y);
    std::vector<double> hermite_x{0.0, 1.0};
    std::vector<double> hermite_y{0.0, 1.0};
    std::vector<double> hermite_dy{0.0, 2.0};
    MatCal::Algorithm::Insert::Hermite hermite(hermite_x, hermite_y, hermite_dy);
    if (lagrange.calculate(1.5) < 3.24 || quotient.calculate(1.5) < 3.24 ||
        finite.calculate(1.5) < 2.24 || hermite.calculate(0.5) < 0.24) {
        return 4;
    }

    auto root = MatCal::Algorithm::Iteration::Bisection::solve(
        [](double x) { return x * x - 4.0; }, 0.0, 3.0);
    if (!(root > 1.999 && root < 2.001)) {
        return 5;
    }

    std::vector<MatCal::Algorithm::Iteration::NewtonForEquations::Function> funcs;
    funcs.push_back([](const std::vector<double>& values) { return values[0] + values[1] - 3.0; });
    funcs.push_back([](const std::vector<double>& values) { return values[0] - values[1] - 1.0; });
    auto system_root = MatCal::Algorithm::Iteration::NewtonForEquations::solve(2, funcs, {0.0, 0.0});
    if (!system_root.converged) {
        return 6;
    }

    auto least_square = MatCal::Algorithm::Basics::Least_Square::solve(1, xs, ys);
    if (least_square.msg != "success!") {
        return 7;
    }

    double theta_out = 0.0;
    double omega_out = 0.0;
    MatCal::Algorithm::Basics::Integrate::RK4::step2(
        [](double theta, double omega, double& dtheta, double& domega) {
            dtheta = omega;
            domega = -theta;
        },
        1.0,
        0.0,
        0.1,
        theta_out,
        omega_out);
    return theta_out > 0.99 && theta_out < 1.0 && omega_out < 0.0 ? 0 : 8;
}
