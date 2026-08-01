#include "MatCal/Interpolation/LinearInterpolator.hpp"
#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/SkylineLdlt.hpp"
#include "MatCal/Linalg/SymmetricSkylineMatrix.hpp"
#include "MatCal/Linalg/Vector.hpp"
#include "MatCal/Polynomial/Polynomial.hpp"
#include "MatCal/Roots/Roots.hpp"

#include <cmath>
#include <iostream>

namespace {

bool close(double a, double b, double tol = 1.0e-10) {
    return std::abs(a - b) <= tol;
}

} // namespace

int main() {
    using namespace MatCal::Linalg;

    DenseMatrix a{{4.0, 1.0}, {1.0, 3.0}};
    Vector b{1.0, 2.0};
    SolverResult direct = solve_dense_partial_pivot(a, b);
    if (!direct.success() || !close(direct.solution[0], 1.0 / 11.0) || !close(direct.solution[1], 7.0 / 11.0)) {
        std::cerr << "dense solve documentation example failed\n";
        return 1;
    }

    SymmetricSkylineMatrix skyline({0, 0});
    skyline.set(0, 0, 4.0);
    skyline.set(1, 0, 1.0);
    skyline.set(1, 1, 3.0);
    auto factorization = factorize_skyline_ldlt(skyline);
    if (!factorization.success()) {
        std::cerr << "skyline factorization documentation example failed\n";
        return 1;
    }
    SolverResult skyline_result = factorization.factorization.solve(b);
    if (!skyline_result.success() || !close(skyline_result.solution[0], direct.solution[0]) ||
        !close(skyline_result.solution[1], direct.solution[1])) {
        std::cerr << "skyline solve documentation example failed\n";
        return 1;
    }

    MatCal::Polynomial::Polynomial p{1.0, 0.0, 1.0};
    if (!close(p.evaluate(3.0), 10.0)) {
        std::cerr << "polynomial documentation example failed\n";
        return 1;
    }

    auto root = MatCal::Roots::solve_bisection([](double x) { return x * x - 2.0; }, 0.0, 2.0);
    if (!root.success() || !close(root.value, std::sqrt(2.0), 1.0e-8)) {
        std::cerr << "roots documentation example failed\n";
        return 1;
    }

    MatCal::Interpolation::LinearInterpolator line({{0.0, 1.0}, {2.0, 5.0}});
    if (!close(line.evaluate(0.5), 2.0)) {
        std::cerr << "interpolation documentation example failed\n";
        return 1;
    }

    return 0;
}
