#include <cmath>
#include <limits>
#include <vector>

#include "MatCal/Linalg/SkylineLdlt.hpp"
#include "MatCal/Linalg/SolverTypes.hpp"
#include "MatCal/Linalg/SymmetricSkylineMatrix.hpp"
#include "MatCal/Linalg/Vector.hpp"

namespace {

bool near(double actual, double expected, double tolerance) {
    return std::abs(actual - expected) <= tolerance;
}

} // namespace

int main() {
    using MatCal::Linalg::SolverOptions;
    using MatCal::Linalg::SymmetricSkylineMatrix;
    using MatCal::Linalg::Vector;

    SymmetricSkylineMatrix stiffness = SymmetricSkylineMatrix::from_profile({0, 0, 1});
    stiffness.add(0, 0, 4.0);
    stiffness.add(1, 0, -1.0);
    stiffness.add(1, 1, 4.0);
    stiffness.add(2, 1, -1.0);
    stiffness.add(2, 2, 3.0);

    Vector expected{1.0, 2.0, 3.0};
    Vector rhs = stiffness.multiply(expected);

    SolverOptions options;
    const double c = 256.0 * std::numeric_limits<double>::epsilon() * static_cast<double>(stiffness.size());
    options.absolute_tolerance = c;
    options.relative_tolerance = c;
    options.pivot_factor = 1.0;

    auto factor = MatCal::Linalg::factorize_skyline_ldlt(stiffness, options);
    if (!factor.success()) {
        return 1;
    }
    if (factor.factorization.storage_size() >= stiffness.size() * stiffness.size()) {
        return 2;
    }
    if (factor.metrics.matrix_scale != stiffness.matrix_scale()) {
        return 3;
    }

    auto solved = factor.factorization.solve(rhs, options);
    if (!solved.success()) {
        return 4;
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (!near(solved.solution[i], expected[i], 1.0e-11)) {
            return 5;
        }
    }

    auto solved_again = factor.factorization.solve(stiffness.multiply(Vector{2.0, -1.0, 0.5}), options);
    if (!solved_again.success()) {
        return 6;
    }
    if (solved_again.metrics.factorization_operation_count != factor.metrics.factorization_operation_count) {
        return 7;
    }
    if (!near(solved.metrics.absolute_residual_norm, 0.0, solved.metrics.residual_acceptance_tolerance)) {
        return 8;
    }
    if (solved.metrics.relative_residual_norm < 0.0 || !std::isfinite(solved.metrics.relative_residual_norm)) {
        return 9;
    }
    return 0;
}
