#include <limits>

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

} // namespace

int main() {
    stationary_solvers_converge();
    scale_invariance();
    failures_are_structured();
    return finish("MatCal::Linalg iterative solver tests");
}
