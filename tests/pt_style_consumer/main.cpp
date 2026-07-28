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

    auto root = MatCal::Algorithm::Iteration::Bisection::solve(
        [](double x) { return x * x - 4.0; }, 0.0, 3.0);
    if (!(root > 1.999 && root < 2.001)) {
        return 4;
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
    return theta_out > 0.99 && theta_out < 1.0 && omega_out < 0.0 ? 0 : 5;
}
