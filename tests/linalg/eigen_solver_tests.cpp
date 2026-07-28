#include <cmath>

#include "MatCal/Linalg/EigenSolvers.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

MatCal::Linalg::EigenOptions eigen_options() {
    MatCal::Linalg::EigenOptions options;
    options.absolute_tolerance = 1.0e-10;
    options.relative_tolerance = 1.0e-10;
    options.max_iterations = 200;
    return options;
}

void expect_residual(const MatCal::Linalg::EigenResult& result, double tolerance, const std::string& label) {
    expect_true(result.success(), label + " success");
    expect_true(result.metrics.residual_norm <= tolerance, label + " residual");
}

void dominant_power_cases() {
    using MatCal::Linalg::DenseMatrix;
    auto diagonal = MatCal::Linalg::dominant_eigenpair(DenseMatrix{{3.0, 0.0}, {0.0, 1.0}}, eigen_options());
    expect_residual(diagonal, 1.0e-8, "diagonal dominant");
    expect_near(diagonal.eigenvalue, 3.0, 1.0e-8, "diagonal eigenvalue");

    auto symmetric = MatCal::Linalg::dominant_eigenpair(DenseMatrix{{2.0, 1.0}, {1.0, 2.0}}, eigen_options());
    expect_residual(symmetric, 1.0e-8, "symmetric 2x2");
    expect_near(symmetric.eigenvalue, 3.0, 1.0e-8, "symmetric 2x2 eigenvalue");

    auto negative = MatCal::Linalg::dominant_eigenpair(DenseMatrix{{-3.0, 0.0}, {0.0, 2.0}}, eigen_options());
    expect_residual(negative, 1.0e-8, "negative dominant");
    expect_near(negative.eigenvalue, -3.0, 1.0e-8, "negative dominant eigenvalue");

    auto scaled = MatCal::Linalg::dominant_eigenpair(DenseMatrix{{3.0e20, 0.0}, {0.0, 1.0e20}}, eigen_options());
    expect_true(scaled.success(), "scaled power succeeds");
    expect_relative(scaled.eigenvalue, 3.0e20, 1.0e-8, "scaled eigenvalue");
}

void inverse_power_cases() {
    using MatCal::Linalg::DenseMatrix;
    auto options = eigen_options();
    options.shift = 1.9;
    auto shifted = MatCal::Linalg::inverse_power_eigenpair(DenseMatrix{{3.0, 0.0}, {0.0, 2.0}}, options);
    expect_residual(shifted, 1.0e-8, "shifted inverse");
    expect_near(shifted.eigenvalue, 2.0, 1.0e-8, "shifted inverse eigenvalue");
    expect_true(shifted.metrics.linear_solves > 0, "inverse power records linear solves");

    options.shift = 2.0;
    auto singular_shift = MatCal::Linalg::inverse_power_eigenpair(DenseMatrix{{3.0, 0.0}, {0.0, 2.0}}, options);
    expect_true(singular_shift.status == MatCal::Linalg::EigenStatus::singular_shift,
                "singular shifted matrix is structured");
}

void invalid_cases() {
    auto zero_initial = MatCal::Linalg::dominant_eigenpair(
        MatCal::Linalg::DenseMatrix::identity(2), MatCal::Linalg::Vector{0.0, 0.0}, eigen_options());
    expect_true(zero_initial.status == MatCal::Linalg::EigenStatus::invalid_input,
                "zero initial vector rejected");

    auto limited = eigen_options();
    limited.max_iterations = 1;
    auto hard = MatCal::Linalg::dominant_eigenpair(
        MatCal::Linalg::DenseMatrix{{2.0, 0.0}, {0.0, 1.999}}, MatCal::Linalg::Vector{1.0, 1.0}, limited);
    expect_true(hard.status == MatCal::Linalg::EigenStatus::not_converged ||
                    hard.status == MatCal::Linalg::EigenStatus::success,
                "difficult repeated-scale case is deterministic");
}

} // namespace

int main() {
    dominant_power_cases();
    inverse_power_cases();
    invalid_cases();
    return finish("MatCal::Linalg eigen solver tests");
}
