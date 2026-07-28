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
    return root > 1.999 && root < 2.001 ? 0 : 4;
}
