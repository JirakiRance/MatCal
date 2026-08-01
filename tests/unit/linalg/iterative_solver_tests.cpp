#include <cmath>
#include <limits>
#include <random>

#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/IterativeSolvers.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

MatCal::Linalg::SolverOptions iterative_options() {
    MatCal::Linalg::SolverOptions options;
    options.absolute_tolerance = 0.0;
    options.relative_tolerance = 1.0e-10;
    options.max_iterations = 1000;
    return options;
}

void expect_solution(const MatCal::Linalg::SolverResult& result, const std::string& label) {
    expect_true(result.success(), label + " success");
    expect_near(result.solution[0], 0.1, 1.0e-7, label + " x0");
    expect_near(result.solution[1], 0.6, 1.0e-7, label + " x1");
    expect_true(result.metrics.absolute_residual_norm < result.metrics.residual_acceptance_tolerance * 1.01,
                label + " residual accepted");
}

void stationary_solvers_converge() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::Vector;
    DenseMatrix a{{4.0, 1.0}, {2.0, 3.0}};
    Vector b{1.0, 2.0};
    auto options = iterative_options();

    expect_solution(MatCal::Linalg::solve_jacobi(a, b, options), "Jacobi");
    auto gs = MatCal::Linalg::solve_gauss_seidel(a, b, options);
    auto sor1 = MatCal::Linalg::solve_sor(a, b, 1.0, options);
    expect_solution(gs, "Gauss-Seidel");
    expect_solution(sor1, "SOR omega=1");
    expect_near(gs.solution[0], sor1.solution[0], 1.0e-10, "SOR omega=1 matches GS x0");
    expect_near(gs.solution[1], sor1.solution[1], 1.0e-10, "SOR omega=1 matches GS x1");

    expect_solution(MatCal::Linalg::solve_sor(a, b, 1.2, options), "SOR omega=1.2");
    expect_solution(MatCal::Linalg::solve_sor(a, b, 1.5, options), "SOR omega=1.5");

    auto direct = MatCal::Linalg::solve_dense_partial_pivot(a, b);
    expect_near(gs.solution[0], direct.solution[0], 1.0e-8, "GS matches direct x0");
    expect_near(gs.solution[1], direct.solution[1], 1.0e-8, "GS matches direct x1");
}

void scale_invariance() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::Vector;
    auto options = iterative_options();
    for (double scale : {1.0e-20, 1.0, 1.0e20}) {
        DenseMatrix a{{4.0 * scale, 1.0 * scale}, {2.0 * scale, 3.0 * scale}};
        Vector b{1.0 * scale, 2.0 * scale};
        auto result = MatCal::Linalg::solve_gauss_seidel(a, b, options);
        expect_true(result.success(), "GS scaled success");
        expect_near(result.solution[0], 0.1, 1.0e-7, "GS scaled x0");
        expect_near(result.solution[1], 0.6, 1.0e-7, "GS scaled x1");
    }
}

void failures_are_structured() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::SolverStatus;
    using MatCal::Linalg::Vector;
    auto options = iterative_options();

    auto zero_diag = MatCal::Linalg::solve_jacobi(DenseMatrix{{0.0, 1.0}, {1.0, 2.0}}, Vector{1.0, 1.0}, options);
    expect_true(zero_diag.status == SolverStatus::singular, "zero diagonal is singular status");

    auto invalid_omega = MatCal::Linalg::solve_sor(DenseMatrix::identity(1), Vector{1.0}, 2.0, options);
    expect_true(invalid_omega.status == SolverStatus::invalid_input, "invalid omega structured");

    options.max_iterations = 2;
    auto divergent = MatCal::Linalg::solve_jacobi(DenseMatrix{{1.0, 3.0}, {3.0, 1.0}}, Vector{1.0, 1.0}, options);
    expect_true(divergent.status == SolverStatus::not_converged, "non-convergence structured");
    expect_true(divergent.solution.empty(), "non-convergence returns no partial solution");

    auto nonfinite = MatCal::Linalg::solve_gauss_seidel(
        DenseMatrix{{1.0, std::numeric_limits<double>::quiet_NaN()}, {0.0, 1.0}},
        Vector{1.0, 1.0}, iterative_options());
    expect_true(nonfinite.status == SolverStatus::non_finite_input, "non-finite input structured");
}

void edge_contracts() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::SolverStatus;
    using MatCal::Linalg::Vector;
    auto options = iterative_options();

    auto empty = MatCal::Linalg::solve_gauss_seidel(DenseMatrix(0, 0), Vector{}, options);
    expect_true(empty.success(), "empty system succeeds as empty solution");
    expect_true(empty.solution.empty(), "empty system returns empty solution");

    auto zero_matrix = MatCal::Linalg::solve_jacobi(DenseMatrix{{0.0}}, Vector{0.0}, options);
    expect_true(zero_matrix.status == SolverStatus::singular, "zero matrix is singular even with zero RHS");
    expect_true(zero_matrix.solution.empty(), "singular zero matrix returns no partial solution");

    auto one_by_one = MatCal::Linalg::solve_gauss_seidel(DenseMatrix{{4.0}}, Vector{8.0}, options);
    expect_true(one_by_one.success(), "1x1 GS success");
    expect_near(one_by_one.solution[0], 2.0, 1.0e-12, "1x1 solution");

    auto zero_rhs = MatCal::Linalg::solve_gauss_seidel(
        DenseMatrix::identity(2), Vector{0.0, 0.0}, Vector{10.0, -4.0}, options);
    expect_true(zero_rhs.success(), "zero RHS converges from nonzero initial guess");
    expect_near(zero_rhs.solution[0], 0.0, 0.0, "zero RHS x0");
    expect_near(zero_rhs.solution[1], 0.0, 0.0, "zero RHS x1");
    expect_near(zero_rhs.metrics.solution_scale, 0.0, 0.0, "zero RHS final solution scale");
    expect_near(zero_rhs.metrics.absolute_residual_norm,
                MatCal::Linalg::residual_norm_inf(DenseMatrix::identity(2), zero_rhs.solution, Vector{0.0, 0.0}),
                0.0,
                "metrics residual matches final solution");

    auto bad_initial_size = MatCal::Linalg::solve_jacobi(
        DenseMatrix::identity(2), Vector{1.0, 1.0}, Vector{0.0}, options);
    expect_true(bad_initial_size.status == SolverStatus::dimension_mismatch, "initial guess size checked");

    auto bad_initial_finite = MatCal::Linalg::solve_jacobi(
        DenseMatrix::identity(1), Vector{1.0}, Vector{std::numeric_limits<double>::infinity()}, options);
    expect_true(bad_initial_finite.status == SolverStatus::non_finite_input, "initial guess finite checked");

    auto invalid_options = options;
    invalid_options.max_iterations = 0;
    auto max_zero = MatCal::Linalg::solve_gauss_seidel(DenseMatrix::identity(1), Vector{1.0}, invalid_options);
    expect_true(max_zero.status == SolverStatus::invalid_input, "max_iterations=0 invalid");

    auto tiny_step = options;
    tiny_step.max_iterations = 1;
    auto nearly_no_update = MatCal::Linalg::solve_sor(DenseMatrix::identity(1), Vector{1.0}, 1.0e-12, tiny_step);
    expect_true(nearly_no_update.status == SolverStatus::not_converged,
                "tiny SOR step does not masquerade as residual success");
    expect_true(nearly_no_update.solution.empty(), "tiny-step failure returns no partial solution");

    auto near_two = MatCal::Linalg::solve_sor(
        DenseMatrix::identity(1), Vector{1.0}, std::nextafter(2.0, 0.0), Vector{0.0}, tiny_step);
    expect_true(near_two.status == SolverStatus::not_converged ||
                    near_two.status == SolverStatus::success,
                "omega just below 2 remains valid and structured");

    const double big = std::numeric_limits<double>::max();
    auto overflow_residual = MatCal::Linalg::solve_jacobi(
        DenseMatrix{{big, big}, {0.0, big}}, Vector{0.0, 0.0}, Vector{1.0, 1.0}, options);
    expect_true(overflow_residual.status == SolverStatus::breakdown,
                "initial residual overflow is breakdown");
    expect_true(overflow_residual.solution.empty(), "residual overflow returns no partial solution");
}

void inputs_are_not_modified_and_random_differential() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::Vector;
    std::mt19937 rng(7131);
    std::uniform_real_distribution<double> dist(-0.25, 0.25);
    auto options = iterative_options();
    options.relative_tolerance = 1.0e-11;
    options.max_iterations = 2000;

    for (int trial = 0; trial < 6; ++trial) {
        DenseMatrix a(4, 4);
        Vector expected{1.0 + trial, -0.5, 0.25, 2.0};
        for (std::size_t r = 0; r < 4; ++r) {
            double row_sum = 0.0;
            for (std::size_t c = 0; c < 4; ++c) {
                if (r == c) {
                    continue;
                }
                const double value = dist(rng);
                a(r, c) = value;
                row_sum += std::abs(value);
            }
            a(r, r) = row_sum + 2.0 + 0.1 * static_cast<double>(trial);
        }
        Vector b = a.multiply(expected);
        DenseMatrix a_before = a;
        Vector b_before = b;
        auto gs = MatCal::Linalg::solve_gauss_seidel(a, b, options);
        auto direct = MatCal::Linalg::solve_dense_partial_pivot(a, b);
        expect_true(gs.success(), "random diagonal-dominant GS success");
        expect_true(direct.success(), "random diagonal-dominant direct success");
        for (std::size_t i = 0; i < 4; ++i) {
            expect_near(gs.solution[i], direct.solution[i], 1.0e-8, "random GS/direct differential");
            expect_near(b[i], b_before[i], 0.0, "iterative RHS unchanged");
            for (std::size_t j = 0; j < 4; ++j) {
                expect_near(a(i, j), a_before(i, j), 0.0, "iterative matrix unchanged");
            }
        }
    }
}

} // namespace

int main() {
    stationary_solvers_converge();
    scale_invariance();
    failures_are_structured();
    edge_contracts();
    inputs_are_not_modified_and_random_differential();
    return finish("MatCal::Linalg iterative solver tests");
}
