#include <cmath>
#include <limits>
#include <vector>

#include "test_support.hpp"
#include "Iteration.hpp"
#include "MatCal/Nonlinear/Nonlinear.hpp"

namespace {

using matcal_test::expect_near;
using matcal_test::expect_true;
using matcal_test::finish;

MatCal::Nonlinear::NonlinearOptions options() {
    MatCal::Nonlinear::NonlinearOptions opts;
    opts.absolute_tolerance = 1.0e-11;
    opts.relative_tolerance = 1.0e-11;
    opts.finite_difference_step = 1.0e-6;
    opts.max_iterations = 25;
    return opts;
}

void linear_system_with_analytic_jacobian() {
    auto residual = [](const std::vector<double>& x) {
        return std::vector<double>{x[0] + x[1] - 3.0, x[0] - x[1] - 1.0};
    };
    auto jacobian = [](const std::vector<double>&) {
        MatCal::Linalg::DenseMatrix j(2, 2);
        j(0, 0) = 1.0; j(0, 1) = 1.0;
        j(1, 0) = 1.0; j(1, 1) = -1.0;
        return j;
    };

    auto result = MatCal::Nonlinear::solve_newton_system(residual, {0.0, 0.0}, jacobian, options());
    expect_true(result.success(), "analytic Newton solves linear system");
    expect_near(result.solution[0], 2.0, 1.0e-12, "analytic Newton x0");
    expect_near(result.solution[1], 1.0, 1.0e-12, "analytic Newton x1");
    expect_true(result.metrics.linear_solves == 1, "analytic Newton records linear solve");
}

void finite_difference_jacobian_matches_analytic_solution() {
    auto residual = [](const std::vector<double>& x) {
        return std::vector<double>{x[0] * x[0] + x[1] * x[1] - 1.0, x[0] - x[1]};
    };
    auto result = MatCal::Nonlinear::solve_newton_system_finite_difference(residual, {0.8, 0.7}, options());
    const double expected = std::sqrt(0.5);
    expect_true(result.success(), "finite difference Newton solves circle-line system");
    expect_near(result.solution[0], expected, 1.0e-8, "finite difference Newton x0");
    expect_near(result.solution[1], expected, 1.0e-8, "finite difference Newton x1");
    expect_true(result.metrics.function_evaluations > result.metrics.iterations, "finite difference counts residual calls");
}

void failure_cases_are_structured() {
    auto residual = [](const std::vector<double>& x) {
        return std::vector<double>{x[0] * x[0], x[1]};
    };
    auto singular_jacobian = [](const std::vector<double>& x) {
        MatCal::Linalg::DenseMatrix j(2, 2);
        j(0, 0) = 2.0 * x[0]; j(0, 1) = 0.0;
        j(1, 0) = 0.0;        j(1, 1) = 1.0;
        return j;
    };
    auto singular = MatCal::Nonlinear::solve_newton_system(residual, {0.0, 1.0}, singular_jacobian, options());
    expect_true(singular.status == MatCal::Nonlinear::NonlinearStatus::singular_jacobian,
                "singular Jacobian is not reported as success");
    expect_true(singular.solution.empty(), "singular Jacobian returns no partial solution");

    auto mismatch = MatCal::Nonlinear::solve_newton_system(
        [](const std::vector<double>& x) { return std::vector<double>{x[0]}; },
        {1.0, 2.0}, singular_jacobian, options());
    expect_true(mismatch.status == MatCal::Nonlinear::NonlinearStatus::dimension_mismatch,
                "residual dimension mismatch is structured");

    auto nonfinite = MatCal::Nonlinear::solve_newton_system(
        [](const std::vector<double>&) { return std::vector<double>{std::numeric_limits<double>::quiet_NaN()}; },
        {1.0},
        [](const std::vector<double>&) {
            MatCal::Linalg::DenseMatrix j(1, 1);
            j(0, 0) = 1.0;
            return j;
        },
        options());
    expect_true(nonfinite.diagnostic.reason == MatCal::Nonlinear::NonlinearReason::non_finite_residual,
                "non-finite residual has machine-readable reason");

    auto limited_options = options();
    limited_options.max_iterations = 1;
    auto no_root = MatCal::Nonlinear::solve_newton_system(
        [](const std::vector<double>& x) { return std::vector<double>{std::exp(x[0]) + 2.0}; },
        {0.0},
        [](const std::vector<double>& x) {
            MatCal::Linalg::DenseMatrix j(1, 1);
            j(0, 0) = std::exp(x[0]);
            return j;
        },
        limited_options);
    expect_true(no_root.status == MatCal::Nonlinear::NonlinearStatus::not_converged,
                "maximum iteration failure is not pseudo-success");
}

void legacy_newton_for_equations_delegates_core() {
    using Legacy = MatCal::Algorithm::Iteration::NewtonForEquations;
    std::vector<Legacy::Function> funcs;
    funcs.push_back([](const std::vector<double>& x) { return x[0] + x[1] - 3.0; });
    funcs.push_back([](const std::vector<double>& x) { return x[0] - x[1] - 1.0; });

    auto legacy = Legacy::solve(2, funcs, {0.0, 0.0}, 20, 1.0e-6);
    auto core = MatCal::Nonlinear::solve_newton_system_finite_difference(
        [](const std::vector<double>& x) {
            return std::vector<double>{x[0] + x[1] - 3.0, x[0] - x[1] - 1.0};
        },
        {0.0, 0.0}, options());

    expect_true(legacy.converged, "legacy NewtonForEquations converges through core");
    expect_true(core.success(), "core comparison converges");
    expect_near(legacy.root[0], core.solution[0], 1.0e-8, "legacy/core Newton x0 differential");
    expect_near(legacy.root[1], core.solution[1], 1.0e-8, "legacy/core Newton x1 differential");
}

} // namespace

int main() {
    linear_system_with_analytic_jacobian();
    finite_difference_jacobian_matches_analytic_solution();
    failure_cases_are_structured();
    legacy_newton_for_equations_delegates_core();
    return finish("MatCal::Nonlinear tests");
}
