#include "MatCal/Linalg/DenseSolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace MatCal::Linalg {
namespace {

double saturated_product(double left, double right) noexcept {
    if (left == 0.0 || right == 0.0) {
        return 0.0;
    }
    double max = std::numeric_limits<double>::max();
    if (left > max / right) {
        return max;
    }
    return left * right;
}

double saturated_sum(double left, double right) noexcept {
    double max = std::numeric_limits<double>::max();
    if (left > max - right) {
        return max;
    }
    return left + right;
}

double max_abs_vector_scale(const Vector& vector) noexcept {
    double scale = 0.0;
    for (double value : vector) {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

double max_abs_matrix_scale(const DenseMatrix& matrix) noexcept {
    double scale = 0.0;
    for (double value : matrix.span()) {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

SolverDiagnostic diagnostic(SolverStatus status,
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

SolverResult make_result(SolverStatus status,
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
    result.method = "dense_partial_pivot";
    result.implementation = "MatCal::Linalg reference";
    result.metrics = metrics;
    if (!message.empty()) {
        result.diagnostics.push_back(diagnostic(status,
                                                std::move(code),
                                                std::move(reason),
                                                std::move(phase),
                                                row,
                                                column,
                                                value,
                                                metrics.matrix_scale,
                                                tolerance,
                                                std::move(message)));
    }
    return result;
}

SolverResult breakdown_result(std::string reason,
                              std::string phase,
                              SolverMetrics metrics,
                              std::size_t row,
                              std::size_t column,
                              double value) {
    return make_result(SolverStatus::breakdown,
                       "non_finite_intermediate",
                       std::move(reason),
                       std::move(phase),
                       "Dense solve produced a non-finite intermediate value",
                       metrics,
                       row,
                       column,
                       value);
}

void swap_rows(std::vector<double>& values,
               DenseMatrix::size_type cols,
               DenseMatrix::size_type row_a,
               DenseMatrix::size_type row_b) {
    if (row_a == row_b) {
        return;
    }
    for (DenseMatrix::size_type c = 0; c < cols; ++c) {
        std::swap(values[row_a * cols + c], values[row_b * cols + c]);
    }
}

} // namespace

double residual_norm_inf(const DenseMatrix& a, const Vector& x, const Vector& b) {
    if (a.cols() != x.size() || a.rows() != b.size()) {
        throw std::invalid_argument("residual dimension mismatch");
    }

    double max_residual = 0.0;
    for (DenseMatrix::size_type r = 0; r < a.rows(); ++r) {
        double sum = 0.0;
        for (DenseMatrix::size_type c = 0; c < a.cols(); ++c) {
            sum += a(r, c) * x[c];
        }
        max_residual = std::max(max_residual, std::abs(sum - b[r]));
    }
    return max_residual;
}

SolverResult solve_dense_partial_pivot(const DenseMatrix& a,
                                       const Vector& b,
                                       const SolverOptions& options) {
    if (!options.valid()) {
        return make_result(SolverStatus::invalid_input,
                           "invalid_options",
                           "invalid_solver_options",
                           "validation",
                           "Invalid solver options");
    }
    if (a.rows() != a.cols()) {
        return make_result(SolverStatus::dimension_mismatch,
                           "matrix_not_square",
                           "dimension_mismatch",
                           "validation",
                           "Dense solve requires a square matrix");
    }
    if (a.rows() != b.size()) {
        return make_result(SolverStatus::dimension_mismatch,
                           "rhs_size_mismatch",
                           "dimension_mismatch",
                           "validation",
                           "Right-hand side size does not match matrix rows");
    }
    if (!a.all_finite() || !b.all_finite()) {
        return make_result(SolverStatus::non_finite_input,
                           "non_finite_input",
                           "input_not_finite",
                           "validation",
                           "Dense solve inputs must be finite");
    }

    const auto n = a.rows();
    SolverResult result;
    result.status = SolverStatus::success;
    result.method = "dense_partial_pivot";
    result.implementation = "MatCal::Linalg reference";
    result.solution = Vector(n);
    result.metrics.matrix_scale = max_abs_matrix_scale(a);
    result.metrics.rhs_scale = max_abs_vector_scale(b);
    result.metrics.pivot_tolerance_used = options.pivot_tolerance(result.metrics.matrix_scale);

    if (n == 0) {
        result.metrics.residual_norm = 0.0;
        result.metrics.absolute_residual_norm = 0.0;
        result.metrics.relative_residual_norm = 0.0;
        result.metrics.residual_acceptance_tolerance = options.comparison_tolerance(0.0);
        result.metrics.minimum_abs_pivot = 0.0;
        return result;
    }

    if (!std::isfinite(result.metrics.matrix_scale) ||
        !std::isfinite(result.metrics.rhs_scale) ||
        !std::isfinite(result.metrics.pivot_tolerance_used)) {
        return breakdown_result("scale_not_finite",
                                "scale",
                                result.metrics,
                                SolverDiagnostic::invalid_index(),
                                SolverDiagnostic::invalid_index(),
                                result.metrics.matrix_scale);
    }

    std::vector<double> work = a.values();
    std::vector<double> rhs = b.values();

    for (DenseMatrix::size_type k = 0; k < n; ++k) {
        DenseMatrix::size_type pivot_row = k;
        double pivot_abs = std::abs(work[k * n + k]);
        for (DenseMatrix::size_type r = k + 1; r < n; ++r) {
            double candidate = std::abs(work[r * n + k]);
            if (!std::isfinite(candidate)) {
                return breakdown_result("pivot_search_not_finite", "pivot_search", result.metrics, r, k, candidate);
            }
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot_row = r;
            }
        }

        if (!std::isfinite(pivot_abs)) {
            return breakdown_result("pivot_not_finite", "pivot_search", result.metrics, pivot_row, k, pivot_abs);
        }
        result.metrics.minimum_abs_pivot = std::min(result.metrics.minimum_abs_pivot, pivot_abs);

        if (pivot_abs <= result.metrics.pivot_tolerance_used) {
            return make_result(SolverStatus::singular,
                               "pivot_too_small",
                               "pivot_below_tolerance",
                               "factorization",
                               "Pivot is zero or below pivot tolerance",
                               result.metrics,
                               pivot_row,
                               k,
                               pivot_abs,
                               result.metrics.pivot_tolerance_used);
        }

        swap_rows(work, n, k, pivot_row);
        std::swap(rhs[k], rhs[pivot_row]);

        double pivot = work[k * n + k];
        if (!std::isfinite(pivot)) {
            return breakdown_result("pivot_not_finite", "factorization", result.metrics, k, k, pivot);
        }
        for (DenseMatrix::size_type r = k + 1; r < n; ++r) {
            double factor = work[r * n + k] / pivot;
            if (!std::isfinite(factor)) {
                return breakdown_result("factor_not_finite", "elimination", result.metrics, r, k, factor);
            }
            work[r * n + k] = 0.0;
            for (DenseMatrix::size_type c = k + 1; c < n; ++c) {
                work[r * n + c] -= factor * work[k * n + c];
                ++result.metrics.operation_count;
                if (!std::isfinite(work[r * n + c])) {
                    return breakdown_result("matrix_update_not_finite", "elimination", result.metrics, r, c, work[r * n + c]);
                }
            }
            rhs[r] -= factor * rhs[k];
            if (!std::isfinite(rhs[r])) {
                return breakdown_result("rhs_update_not_finite", "elimination", result.metrics, r, SolverDiagnostic::invalid_index(), rhs[r]);
            }
        }
        result.metrics.iterations = k + 1;
    }
    result.metrics.factorization_operation_count = result.metrics.operation_count;

    for (DenseMatrix::size_type i = n; i-- > 0;) {
        double sum = rhs[i];
        for (DenseMatrix::size_type c = i + 1; c < n; ++c) {
            sum -= work[i * n + c] * result.solution[c];
            ++result.metrics.operation_count;
            if (!std::isfinite(sum)) {
                return breakdown_result("back_substitution_sum_not_finite", "back_substitution", result.metrics, i, c, sum);
            }
        }
        double diag = work[i * n + i];
        double diag_abs = std::abs(diag);
        if (!std::isfinite(diag_abs)) {
            return breakdown_result("back_substitution_pivot_not_finite", "back_substitution", result.metrics, i, i, diag);
        }
        result.metrics.minimum_abs_pivot = std::min(result.metrics.minimum_abs_pivot, diag_abs);
        if (diag_abs <= result.metrics.pivot_tolerance_used) {
            return make_result(SolverStatus::singular,
                               "pivot_too_small",
                               "back_substitution_pivot_below_tolerance",
                               "back_substitution",
                               "Back substitution pivot is below tolerance",
                               result.metrics,
                               i,
                               i,
                               diag_abs,
                               result.metrics.pivot_tolerance_used);
        }
        result.solution[i] = sum / diag;
        if (!std::isfinite(result.solution[i])) {
            return breakdown_result("solution_not_finite", "back_substitution", result.metrics, i, SolverDiagnostic::invalid_index(), result.solution[i]);
        }
    }

    result.metrics.solution_scale = max_abs_vector_scale(result.solution);
    if (!std::isfinite(result.metrics.solution_scale)) {
        return breakdown_result("solution_scale_not_finite",
                                "solution",
                                result.metrics,
                                SolverDiagnostic::invalid_index(),
                                SolverDiagnostic::invalid_index(),
                                result.metrics.solution_scale);
    }

    double max_residual = 0.0;
    for (DenseMatrix::size_type r = 0; r < a.rows(); ++r) {
        double sum = 0.0;
        for (DenseMatrix::size_type c = 0; c < a.cols(); ++c) {
            sum += a(r, c) * result.solution[c];
            ++result.metrics.operation_count;
            if (!std::isfinite(sum)) {
                return breakdown_result("residual_sum_not_finite", "residual", result.metrics, r, c, sum);
            }
        }
        double residual = sum - b[r];
        if (!std::isfinite(residual)) {
            return breakdown_result("residual_not_finite", "residual", result.metrics, r, SolverDiagnostic::invalid_index(), residual);
        }
        max_residual = std::max(max_residual, std::abs(residual));
    }

    result.metrics.absolute_residual_norm = max_residual;
    result.metrics.residual_norm = max_residual;

    double denominator = saturated_sum(saturated_product(result.metrics.matrix_scale, result.metrics.solution_scale),
                                      result.metrics.rhs_scale);
    if (denominator == 0.0) {
        result.metrics.relative_residual_norm = max_residual == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    } else {
        result.metrics.relative_residual_norm = max_residual / denominator;
    }
    if (!std::isfinite(result.metrics.relative_residual_norm)) {
        return breakdown_result("relative_residual_not_finite", "residual", result.metrics, SolverDiagnostic::invalid_index(), SolverDiagnostic::invalid_index(), result.metrics.relative_residual_norm);
    }

    result.metrics.residual_acceptance_tolerance = options.comparison_tolerance(denominator);
    if (!std::isfinite(result.metrics.residual_acceptance_tolerance)) {
        return breakdown_result("residual_tolerance_not_finite", "residual", result.metrics, SolverDiagnostic::invalid_index(), SolverDiagnostic::invalid_index(), result.metrics.residual_acceptance_tolerance);
    }
    if (result.metrics.absolute_residual_norm > result.metrics.residual_acceptance_tolerance) {
        result.metrics.solve_operation_count =
            result.metrics.operation_count - result.metrics.factorization_operation_count;
        return make_result(SolverStatus::not_converged,
                           "residual_too_large",
                           "residual_above_tolerance",
                           "residual",
                           "Dense solve residual is above acceptance tolerance",
                           result.metrics,
                           SolverDiagnostic::invalid_index(),
                           SolverDiagnostic::invalid_index(),
                           result.metrics.absolute_residual_norm,
                           result.metrics.residual_acceptance_tolerance);
    }

    result.metrics.solve_operation_count =
        result.metrics.operation_count - result.metrics.factorization_operation_count;
    return result;
}

} // namespace MatCal::Linalg
