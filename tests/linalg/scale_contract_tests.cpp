#include <limits>
#include <string>

#include "MatCal/Linalg/DenseSolver.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

MatCal::Linalg::DenseMatrix scaled_matrix(double scale) {
    return MatCal::Linalg::DenseMatrix{{4.0 * scale, 1.0 * scale},
                                       {2.0 * scale, 3.0 * scale}};
}

MatCal::Linalg::Vector scaled_rhs(double scale) {
    return MatCal::Linalg::Vector{1.0 * scale, 2.0 * scale};
}

void expect_solution_01_06(const MatCal::Linalg::SolverResult& result, const std::string& label) {
    expect_true(result.success(), label + " success");
    expect_near(result.solution[0], 0.1, 1e-12, label + " x0");
    expect_near(result.solution[1], 0.6, 1e-12, label + " x1");
    expect_true(result.metrics.absolute_residual_norm <= result.metrics.residual_acceptance_tolerance,
                label + " residual accepted");
    expect_true(result.metrics.matrix_scale > 0.0, label + " matrix scale recorded");
    expect_true(result.metrics.pivot_tolerance_used < result.metrics.minimum_abs_pivot,
                label + " pivot tolerance below minimum pivot");
}

} // namespace

int main() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::SolverOptions;
    using MatCal::Linalg::SolverStatus;
    using MatCal::Linalg::Vector;
    using MatCal::Linalg::residual_norm_inf;
    using MatCal::Linalg::solve_dense_partial_pivot;

    for (double scale : {1e-20, 1.0, 1e20}) {
        auto result = solve_dense_partial_pivot(scaled_matrix(scale), scaled_rhs(scale));
        expect_solution_01_06(result, "scale " + std::to_string(scale));
    }

    DenseMatrix tiny_identity{{1e-20, 0.0}, {0.0, 1e-20}};
    auto tiny_identity_result = solve_dense_partial_pivot(tiny_identity, Vector{1e-20, 2e-20});
    expect_true(tiny_identity_result.success(), "1e-20 identity is nonsingular");
    expect_near(tiny_identity_result.solution[0], 1.0, 1e-12, "tiny identity x0");
    expect_near(tiny_identity_result.solution[1], 2.0, 1e-12, "tiny identity x1");

    DenseMatrix huge_identity{{1e20, 0.0}, {0.0, 1e20}};
    auto huge_identity_result = solve_dense_partial_pivot(huge_identity, Vector{1e20, 2e20});
    expect_true(huge_identity_result.success(), "1e20 identity is nonsingular");
    expect_near(huge_identity_result.solution[0], 1.0, 1e-12, "huge identity x0");
    expect_near(huge_identity_result.solution[1], 2.0, 1e-12, "huge identity x1");

    for (double scale : {1e-20, 1.0, 1e20}) {
        DenseMatrix singular{{scale, 2.0 * scale}, {2.0 * scale, 4.0 * scale}};
        auto result = solve_dense_partial_pivot(singular, Vector{scale, 2.0 * scale});
        expect_true(result.status == SolverStatus::singular, "scaled singular rejected " + std::to_string(scale));
        expect_true(result.solution.empty(), "singular result has no partial solution");
        expect_true(!result.diagnostics.empty() && result.diagnostics[0].code == "pivot_too_small",
                    "singular diagnostic code stable");
    }

    SolverOptions default_options;
    DenseMatrix nearly_dependent{{1.0, 1.0}, {1.0, 1.0 + 1e-14}};
    auto default_near = solve_dense_partial_pivot(nearly_dependent, Vector{2.0, 2.0 + 1e-14}, default_options);
    expect_true(default_near.status == SolverStatus::singular, "default options reject relative near singular");

    SolverOptions strict_zero_abs;
    strict_zero_abs.absolute_tolerance = 0.0;
    strict_zero_abs.relative_tolerance = 0.0;
    auto strict_near = solve_dense_partial_pivot(nearly_dependent, Vector{2.0, 2.0 + 1e-14}, strict_zero_abs);
    expect_true(strict_near.success(), "zero absolute and zero relative tolerance accept nonzero pivot");

    SolverOptions loose;
    loose.absolute_tolerance = 0.0;
    loose.relative_tolerance = 1e-3;
    auto loose_near = solve_dense_partial_pivot(nearly_dependent, Vector{2.0, 2.0 + 1e-14}, loose);
    expect_true(loose_near.status == SolverStatus::singular, "relative tolerance controls near-singular rejection");

    auto zero_zero = solve_dense_partial_pivot(DenseMatrix(2, 2), Vector(2));
    expect_true(zero_zero.status == SolverStatus::singular, "zero matrix zero rhs is singular");
    auto zero_nonzero = solve_dense_partial_pivot(DenseMatrix(2, 2), Vector{1.0, 0.0});
    expect_true(zero_nonzero.status == SolverStatus::singular, "zero matrix nonzero rhs is singular");

    DenseMatrix overflow_case{{1e308, 1e308}, {1e308, -1e308}};
    auto overflow_result = solve_dense_partial_pivot(overflow_case, Vector{1e308, 1.0});
    expect_true(overflow_result.status == SolverStatus::breakdown, "finite input overflow returns breakdown");
    expect_true(overflow_result.solution.empty(), "breakdown result has no partial solution");
    expect_true(!overflow_result.diagnostics.empty() && overflow_result.diagnostics[0].phase == "elimination",
                "breakdown diagnostic phase");

    auto base = solve_dense_partial_pivot(scaled_matrix(1.0), scaled_rhs(1.0));
    double direct_abs = residual_norm_inf(scaled_matrix(1.0), base.solution, scaled_rhs(1.0));
    expect_near(base.metrics.absolute_residual_norm, direct_abs, 1e-15, "absolute residual matches direct computation");
    double denominator = base.metrics.matrix_scale * base.metrics.solution_scale + base.metrics.rhs_scale;
    expect_near(base.metrics.relative_residual_norm,
                denominator == 0.0 ? 0.0 : direct_abs / denominator,
                1e-15,
                "relative residual definition");

    SolverOptions nan_options;
    nan_options.relative_tolerance = std::numeric_limits<double>::quiet_NaN();
    expect_true(!nan_options.valid(), "NaN option invalid");
    SolverOptions inf_options;
    inf_options.absolute_tolerance = std::numeric_limits<double>::infinity();
    expect_true(!inf_options.valid(), "Inf option invalid");
    SolverOptions negative_pivot;
    negative_pivot.pivot_factor = -1.0;
    expect_true(!negative_pivot.valid(), "negative pivot factor invalid");
    SolverOptions zero_iterations;
    zero_iterations.max_iterations = 0;
    expect_true(!zero_iterations.valid(), "zero max iterations invalid");
    SolverOptions zero_pivot_factor;
    zero_pivot_factor.pivot_factor = 0.0;
    expect_true(zero_pivot_factor.valid(), "zero pivot factor is valid exact-pivot mode");

    auto first = solve_dense_partial_pivot(DenseMatrix{{1.0, 2.0}, {2.0, 4.0}}, Vector{1.0, 2.0});
    auto second = solve_dense_partial_pivot(DenseMatrix{{1.0, 2.0}, {2.0, 4.0}}, Vector{1.0, 2.0});
    expect_true(first.status == second.status, "deterministic status");
    expect_true(first.diagnostics[0].code == second.diagnostics[0].code, "deterministic diagnostic code");
    expect_true(first.diagnostics[0].phase == second.diagnostics[0].phase, "deterministic diagnostic phase");

    return finish("MatCal::Linalg scale contract tests");
}
