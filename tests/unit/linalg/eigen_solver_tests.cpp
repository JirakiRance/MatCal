#include <cmath>
#include <limits>
#include <random>

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

    auto tiny_options = eigen_options();
    tiny_options.absolute_tolerance = 0.0;
    auto tiny = MatCal::Linalg::dominant_eigenpair(DenseMatrix{{3.0e-200, 0.0}, {0.0, 1.0e-200}}, tiny_options);
    expect_true(tiny.success(), "tiny scaled power succeeds");
    expect_near(tiny.eigenvalue, 3.0e-200, 1.0e-208, "tiny scaled eigenvalue");
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
    expect_true(singular_shift.eigenvector.empty(), "singular shifted solve returns no eigenvector");

    options.shift = 2.0 + 1.0e-8;
    auto near_shift = MatCal::Linalg::inverse_power_eigenpair(DenseMatrix{{3.0, 0.0}, {0.0, 2.0}}, options);
    expect_residual(near_shift, 1.0e-8, "near shifted inverse");
    expect_near(near_shift.eigenvalue, 2.0, 1.0e-8, "near shifted inverse eigenvalue");
}

void invalid_cases() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::EigenStatus;
    using MatCal::Linalg::Vector;

    auto zero_initial = MatCal::Linalg::dominant_eigenpair(
        DenseMatrix::identity(2), Vector{0.0, 0.0}, eigen_options());
    expect_true(zero_initial.status == EigenStatus::invalid_input, "zero initial vector rejected");
    expect_true(zero_initial.eigenvector.empty(), "zero initial vector returns no eigenvector");

    auto near_zero_initial = MatCal::Linalg::dominant_eigenpair(
        DenseMatrix{{2.0, 0.0}, {0.0, 1.0}},
        Vector{std::numeric_limits<double>::denorm_min(), 0.0},
        eigen_options());
    expect_true(near_zero_initial.success(), "near-zero but nonzero initial vector is normalized");

    auto limited = eigen_options();
    limited.max_iterations = 1;
    auto hard = MatCal::Linalg::dominant_eigenpair(
        DenseMatrix{{2.0, 0.0}, {0.0, 1.999}}, Vector{1.0, 1.0}, limited);
    expect_true(hard.status == EigenStatus::not_converged || hard.status == EigenStatus::success,
                "difficult close-scale case is deterministic");

    auto repeated = eigen_options();
    repeated.absolute_tolerance = 0.0;
    repeated.max_iterations = 4;
    auto equal_magnitude = MatCal::Linalg::dominant_eigenpair(
        DenseMatrix{{2.0, 0.0}, {0.0, -2.0}}, Vector{1.0, 1.0}, repeated);
    expect_true(equal_magnitude.status == EigenStatus::not_converged,
                "repeated dominant magnitude does not report pseudo success");
    expect_true(equal_magnitude.eigenvector.empty(), "not_converged eigen solve returns no eigenvector");

    auto invalid_options = eigen_options();
    invalid_options.max_iterations = 0;
    auto max_zero = MatCal::Linalg::dominant_eigenpair(DenseMatrix::identity(1), invalid_options);
    expect_true(max_zero.status == EigenStatus::invalid_input, "eigen max_iterations=0 invalid");

    auto nonfinite = MatCal::Linalg::dominant_eigenpair(
        DenseMatrix{{1.0, std::numeric_limits<double>::infinity()}, {0.0, 1.0}}, eigen_options());
    expect_true(nonfinite.status == EigenStatus::non_finite_input, "non-finite matrix rejected");
}

void input_and_random_residual_contracts() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::Vector;

    DenseMatrix a{{2.0, 1.0}, {1.0, 2.0}};
    DenseMatrix before = a;
    auto result = MatCal::Linalg::dominant_eigenpair(a, eigen_options());
    expect_true(result.success(), "input preservation power setup success");
    for (std::size_t r = 0; r < a.rows(); ++r) {
        for (std::size_t c = 0; c < a.cols(); ++c) {
            expect_near(a(r, c), before(r, c), 0.0, "power input matrix unchanged");
        }
    }

    std::mt19937 rng(20260728);
    std::uniform_real_distribution<double> dist(-0.05, 0.05);
    for (int trial = 0; trial < 5; ++trial) {
        DenseMatrix sym(3, 3);
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = r; c < 3; ++c) {
                const double value = dist(rng);
                sym(r, c) = value;
                sym(c, r) = value;
            }
        }
        sym(0, 0) += 1.0;
        sym(1, 1) += 2.0;
        sym(2, 2) += 4.0 + static_cast<double>(trial) * 0.1;
        auto eigen = MatCal::Linalg::dominant_eigenpair(sym, eigen_options());
        expect_true(eigen.success(), "fixed-seed symmetric dominant success");
        expect_true(eigen.metrics.residual_norm <= eigen.metrics.residual_acceptance_tolerance * 1.01,
                    "fixed-seed symmetric residual accepted");
    }
}

} // namespace

int main() {
    dominant_power_cases();
    inverse_power_cases();
    invalid_cases();
    input_and_random_residual_contracts();
    return finish("MatCal::Linalg eigen solver tests");
}
