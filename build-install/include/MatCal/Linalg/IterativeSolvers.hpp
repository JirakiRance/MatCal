#ifndef MATCAL_LINALG_ITERATIVE_SOLVERS_HPP
#define MATCAL_LINALG_ITERATIVE_SOLVERS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>

#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/SolverTypes.hpp"
#include "MatCal/Linalg/Vector.hpp"

namespace MatCal::Linalg {

namespace iterative_detail {

inline double saturated_product(double left, double right) noexcept {
    if (left == 0.0 || right == 0.0) {
        return 0.0;
    }
    const double max = std::numeric_limits<double>::max();
    if (std::abs(left) > max / std::abs(right)) {
        return max;
    }
    return left * right;
}

inline double saturated_sum(double left, double right) noexcept {
    const double max = std::numeric_limits<double>::max();
    if (left > max - right) {
        return max;
    }
    return left + right;
}

inline double matrix_scale(const DenseMatrix& matrix) noexcept {
    double scale = 0.0;
    for (double value : matrix.span()) {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

inline double vector_scale(const Vector& vector) noexcept {
    double scale = 0.0;
    for (double value : vector) {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

inline SolverDiagnostic diagnostic(SolverStatus status,
                                   std::string code,
                                   std::string reason,
                                   std::string phase,
                                   std::size_t row,
                                   std::size_t column,
                                   double value,
                                   double scale,
                                   double tolerance,
                                   std::string message) {
    SolverDiagnostic diag;
    diag.status = status;
    diag.code = std::move(code);
    diag.reason = std::move(reason);
    diag.phase = std::move(phase);
    diag.row = row;
    diag.column = column;
    diag.value = value;
    diag.scale = scale;
    diag.tolerance = tolerance;
    diag.message = std::move(message);
    return diag;
}

inline SolverResult make_result(SolverStatus status,
                                std::string method,
                                std::string code,
                                std::string reason,
                                std::string phase,
                                std::string message,
                                SolverMetrics metrics = {},
                                std::size_t row = SolverDiagnostic::invalid_index(),
                                std::size_t column = SolverDiagnostic::invalid_index(),
                                double value = 0.0,
                                double tolerance = 0.0) {
    SolverResult result;
    result.status = status;
    result.method = std::move(method);
    result.implementation = "MatCal::Linalg stationary iteration";
    result.metrics = metrics;
    if (!message.empty()) {
        result.diagnostics.push_back(diagnostic(status, std::move(code), std::move(reason),
                                                std::move(phase), row, column, value,
                                                metrics.matrix_scale, tolerance,
                                                std::move(message)));
    }
    return result;
}

inline void update_residual_metrics(const DenseMatrix& a,
                                    const Vector& x,
                                    const Vector& b,
                                    SolverMetrics& metrics,
                                    const SolverOptions& options) {
    metrics.solution_scale = vector_scale(x);
    metrics.absolute_residual_norm = residual_norm_inf(a, x, b);
    metrics.residual_norm = metrics.absolute_residual_norm;
    const double denominator = saturated_sum(saturated_product(metrics.matrix_scale, metrics.solution_scale),
                                             metrics.rhs_scale);
    metrics.relative_residual_norm =
        denominator == 0.0
            ? (metrics.absolute_residual_norm == 0.0 ? 0.0 : std::numeric_limits<double>::infinity())
            : metrics.absolute_residual_norm / denominator;
    metrics.residual_acceptance_tolerance = options.comparison_tolerance(denominator);
}

enum class StationaryMethod {
    jacobi,
    gauss_seidel,
    sor
};

inline const char* method_name(StationaryMethod method) noexcept {
    switch (method) {
    case StationaryMethod::jacobi: return "jacobi";
    case StationaryMethod::gauss_seidel: return "gauss_seidel";
    case StationaryMethod::sor: return "sor";
    }
    return "unknown";
}

inline SolverResult solve_stationary(const DenseMatrix& a,
                                     const Vector& b,
                                     StationaryMethod method,
                                     double omega,
                                     const SolverOptions& options,
                                     const Vector* initial_guess) {
    const std::string method_string = method_name(method);
    if (!options.valid()) {
        return make_result(SolverStatus::invalid_input, method_string, "invalid_options",
                           "invalid_solver_options", "validation", "Invalid solver options");
    }
    if (method == StationaryMethod::sor && (!std::isfinite(omega) || omega <= 0.0 || omega >= 2.0)) {
        return make_result(SolverStatus::invalid_input, method_string, "invalid_omega",
                           "omega_out_of_range", "validation", "SOR omega must be finite and in (0, 2)");
    }
    if (a.rows() != a.cols()) {
        return make_result(SolverStatus::dimension_mismatch, method_string, "matrix_not_square",
                           "dimension_mismatch", "validation", "Stationary iteration requires a square matrix");
    }
    if (a.rows() != b.size()) {
        return make_result(SolverStatus::dimension_mismatch, method_string, "rhs_size_mismatch",
                           "dimension_mismatch", "validation", "Right-hand side size does not match matrix rows");
    }
    if (initial_guess != nullptr && initial_guess->size() != b.size()) {
        return make_result(SolverStatus::dimension_mismatch, method_string, "initial_guess_size_mismatch",
                           "dimension_mismatch", "validation", "Initial guess size does not match right-hand side");
    }
    if (!a.all_finite() || !b.all_finite() || (initial_guess != nullptr && !initial_guess->all_finite())) {
        return make_result(SolverStatus::non_finite_input, method_string, "non_finite_input",
                           "input_not_finite", "validation", "Stationary iteration inputs must be finite");
    }

    SolverMetrics metrics;
    metrics.matrix_scale = matrix_scale(a);
    metrics.rhs_scale = vector_scale(b);
    metrics.pivot_tolerance_used = options.pivot_tolerance(metrics.matrix_scale);
    metrics.minimum_abs_pivot = std::numeric_limits<double>::infinity();

    const auto n = a.rows();
    Vector x = initial_guess == nullptr ? Vector(n) : *initial_guess;
    if (n == 0) {
        SolverResult result = make_result(SolverStatus::success, method_string, "", "", "", "", metrics);
        result.solution = x;
        return result;
    }

    for (DenseMatrix::size_type i = 0; i < n; ++i) {
        const double diag_abs = std::abs(a(i, i));
        metrics.minimum_abs_pivot = std::min(metrics.minimum_abs_pivot, diag_abs);
        if (!std::isfinite(diag_abs)) {
            return make_result(SolverStatus::breakdown, method_string, "diagonal_not_finite",
                               "non_finite_diagonal", "validation", "Matrix diagonal is non-finite",
                               metrics, i, i, a(i, i));
        }
        if (diag_abs <= metrics.pivot_tolerance_used) {
            return make_result(SolverStatus::singular, method_string, "diagonal_too_small",
                               "diagonal_below_tolerance", "validation",
                               "Matrix diagonal is zero or below tolerance",
                               metrics, i, i, diag_abs, metrics.pivot_tolerance_used);
        }
    }

    update_residual_metrics(a, x, b, metrics, options);
    if (!std::isfinite(metrics.absolute_residual_norm) ||
        !std::isfinite(metrics.relative_residual_norm) ||
        !std::isfinite(metrics.residual_acceptance_tolerance)) {
        return make_result(SolverStatus::breakdown, method_string, "initial_residual_not_finite",
                           "non_finite_residual", "residual",
                           "Stationary iteration initial residual is non-finite",
                           metrics, SolverDiagnostic::invalid_index(), SolverDiagnostic::invalid_index(),
                           metrics.absolute_residual_norm);
    }
    if (metrics.absolute_residual_norm <= metrics.residual_acceptance_tolerance) {
        SolverResult result = make_result(SolverStatus::success, method_string, "", "", "", "", metrics);
        result.solution = x;
        return result;
    }

    for (std::size_t iteration = 1; iteration <= options.max_iterations; ++iteration) {
        Vector x_new = method == StationaryMethod::jacobi ? Vector(n) : x;
        for (DenseMatrix::size_type i = 0; i < n; ++i) {
            double sum = 0.0;
            for (DenseMatrix::size_type j = 0; j < n; ++j) {
                if (j == i) {
                    continue;
                }
                const double source = (method != StationaryMethod::jacobi && j < i) ? x_new[j] : x[j];
                sum += a(i, j) * source;
                ++metrics.operation_count;
                if (!std::isfinite(sum)) {
                    return make_result(SolverStatus::breakdown, method_string, "iteration_sum_not_finite",
                                       "non_finite_intermediate", "iteration",
                                       "Stationary iteration produced a non-finite row sum",
                                       metrics, i, j, sum);
                }
            }
            const double gauss_value = (b[i] - sum) / a(i, i);
            const double next_value =
                method == StationaryMethod::sor ? omega * gauss_value + (1.0 - omega) * x[i] : gauss_value;
            if (!std::isfinite(next_value)) {
                return make_result(SolverStatus::breakdown, method_string, "iterate_not_finite",
                                   "non_finite_iterate", "iteration",
                                   "Stationary iteration produced a non-finite iterate",
                                   metrics, i, SolverDiagnostic::invalid_index(), next_value);
            }
            x_new[i] = next_value;
        }

        x = std::move(x_new);
        metrics.iterations = iteration;
        update_residual_metrics(a, x, b, metrics, options);
        if (!std::isfinite(metrics.absolute_residual_norm) ||
            !std::isfinite(metrics.relative_residual_norm) ||
            !std::isfinite(metrics.residual_acceptance_tolerance)) {
            return make_result(SolverStatus::breakdown, method_string, "residual_not_finite",
                               "non_finite_residual", "residual",
                               "Stationary iteration residual is non-finite",
                               metrics, SolverDiagnostic::invalid_index(), SolverDiagnostic::invalid_index(),
                               metrics.absolute_residual_norm);
        }
        if (metrics.absolute_residual_norm <= metrics.residual_acceptance_tolerance) {
            SolverResult result = make_result(SolverStatus::success, method_string, "", "", "", "", metrics);
            result.solution = x;
            return result;
        }
    }

    return make_result(SolverStatus::not_converged, method_string, "maximum_iterations",
                       "maximum_iterations", "convergence",
                       "Stationary iteration reached maximum iterations", metrics);
}

} // namespace iterative_detail

inline SolverResult solve_jacobi(const DenseMatrix& a,
                                 const Vector& b,
                                 const SolverOptions& options = {}) {
    return iterative_detail::solve_stationary(a, b, iterative_detail::StationaryMethod::jacobi, 1.0, options, nullptr);
}

inline SolverResult solve_jacobi(const DenseMatrix& a,
                                 const Vector& b,
                                 const Vector& initial_guess,
                                 const SolverOptions& options = {}) {
    return iterative_detail::solve_stationary(a, b, iterative_detail::StationaryMethod::jacobi, 1.0, options, &initial_guess);
}

inline SolverResult solve_gauss_seidel(const DenseMatrix& a,
                                       const Vector& b,
                                       const SolverOptions& options = {}) {
    return iterative_detail::solve_stationary(a, b, iterative_detail::StationaryMethod::gauss_seidel, 1.0, options, nullptr);
}

inline SolverResult solve_gauss_seidel(const DenseMatrix& a,
                                       const Vector& b,
                                       const Vector& initial_guess,
                                       const SolverOptions& options = {}) {
    return iterative_detail::solve_stationary(a, b, iterative_detail::StationaryMethod::gauss_seidel, 1.0, options, &initial_guess);
}

inline SolverResult solve_sor(const DenseMatrix& a,
                              const Vector& b,
                              double omega,
                              const SolverOptions& options = {}) {
    return iterative_detail::solve_stationary(a, b, iterative_detail::StationaryMethod::sor, omega, options, nullptr);
}

inline SolverResult solve_sor(const DenseMatrix& a,
                              const Vector& b,
                              double omega,
                              const Vector& initial_guess,
                              const SolverOptions& options = {}) {
    return iterative_detail::solve_stationary(a, b, iterative_detail::StationaryMethod::sor, omega, options, &initial_guess);
}

} // namespace MatCal::Linalg

#endif
