#include <limits>

#include "MatCal/Linalg/DenseSolver.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

void expect_solution_01_06(const MatCal::Linalg::SolverResult& result, const std::string& label) {
    expect_true(result.success(), label + " success");
    expect_near(result.solution[0], 0.1, 1e-12, label + " x0");
    expect_near(result.solution[1], 0.6, 1e-12, label + " x1");
    expect_true(result.metrics.residual_norm < 1e-12, label + " residual");
}

} // namespace

int main() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::SolverOptions;
    using MatCal::Linalg::SolverStatus;
    using MatCal::Linalg::Vector;
    using MatCal::Linalg::residual_norm_inf;
    using MatCal::Linalg::solve_dense_partial_pivot;

    SolverOptions defaults;
    expect_true(defaults.valid(), "default SolverOptions valid");
    expect_near(defaults.comparison_tolerance(10.0), 1.0e-11, 1e-20, "relative default tolerance");

    SolverOptions invalid = defaults;
    invalid.absolute_tolerance = -1.0;
    expect_true(!invalid.valid(), "negative tolerance invalid");
    auto invalid_result = solve_dense_partial_pivot(DenseMatrix::identity(1), Vector{1.0}, invalid);
    expect_true(invalid_result.status == SolverStatus::invalid_input, "invalid options return structured status");

    DenseMatrix a{{4.0, 1.0}, {2.0, 3.0}};
    Vector b{1.0, 2.0};
    auto result = solve_dense_partial_pivot(a, b);
    expect_solution_01_06(result, "dense solve");
    expect_near(a(0, 0), 4.0, 1e-12, "solver does not modify input matrix");
    expect_near(b[0], 1.0, 1e-12, "solver does not modify input rhs");
    expect_near(residual_norm_inf(a, result.solution, b), result.metrics.absolute_residual_norm, 1e-14, "residual helper");
    expect_near(result.metrics.residual_norm, result.metrics.absolute_residual_norm, 0.0, "legacy residual alias");

    auto repeat = solve_dense_partial_pivot(a, b);
    expect_near(repeat.solution[0], result.solution[0], 0.0, "deterministic solution x0");
    expect_near(repeat.solution[1], result.solution[1], 0.0, "deterministic solution x1");

    DenseMatrix row_swap{{0.0, 1.0}, {1.0, 1.0}};
    Vector row_swap_b{2.0, 5.0};
    auto row_swap_result = solve_dense_partial_pivot(row_swap, row_swap_b);
    expect_true(row_swap_result.success(), "row swap solve succeeds");
    expect_near(row_swap_result.solution[0], 3.0, 1e-12, "row swap solution x0");
    expect_near(row_swap_result.solution[1], 2.0, 1e-12, "row swap solution x1");

    auto non_square = solve_dense_partial_pivot(DenseMatrix(2, 3), Vector(2));
    expect_true(non_square.status == SolverStatus::dimension_mismatch, "non-square matrix status");
    auto rhs_mismatch = solve_dense_partial_pivot(DenseMatrix::identity(2), Vector(3));
    expect_true(rhs_mismatch.status == SolverStatus::dimension_mismatch, "rhs mismatch status");

    DenseMatrix singular{{1.0, 2.0}, {2.0, 4.0}};
    auto singular_result = solve_dense_partial_pivot(singular, Vector{1.0, 2.0});
    expect_true(singular_result.status == SolverStatus::singular, "singular matrix status");
    expect_true(!singular_result.diagnostics.empty(), "singular diagnostic exists");
    expect_true(singular_result.diagnostics[0].code == "pivot_too_small", "singular diagnostic code");

    SolverOptions strict = defaults;
    strict.absolute_tolerance = 1e-12;
    strict.relative_tolerance = 0.0;
    DenseMatrix near_singular{{1e-14, 0.0}, {0.0, 1.0}};
    auto near_singular_result = solve_dense_partial_pivot(near_singular, Vector{1.0, 1.0}, strict);
    expect_true(near_singular_result.status == SolverStatus::singular, "near-singular pivot tolerance status");

    DenseMatrix non_finite{{1.0, std::numeric_limits<double>::infinity()}, {0.0, 1.0}};
    auto non_finite_result = solve_dense_partial_pivot(non_finite, Vector{1.0, 1.0});
    expect_true(non_finite_result.status == SolverStatus::non_finite_input, "non-finite matrix status");
    auto non_finite_rhs = solve_dense_partial_pivot(DenseMatrix::identity(2), Vector{1.0, std::numeric_limits<double>::quiet_NaN()});
    expect_true(non_finite_rhs.status == SolverStatus::non_finite_input, "non-finite rhs status");

    auto empty_result = solve_dense_partial_pivot(DenseMatrix(), Vector());
    expect_true(empty_result.success() && empty_result.solution.empty(), "empty 0x0 solve succeeds");

    expect_throw([&] { (void)residual_norm_inf(DenseMatrix::identity(2), Vector(1), Vector(2)); }, "residual dimension mismatch throws");

    return finish("MatCal::Linalg solver tests");
}
