#include "MatCal/Linalg/DenseSolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace MatCal::Linalg {
namespace {

double saturated_product(double left, double right) noexcept {
    if (left == 0.0 || right == 0.0) {
        return 0.0;
    }
    const double max = std::numeric_limits<double>::max();
    if (left > max / right) {
        return max;
    }
    return left * right;
}

double saturated_sum(double left, double right) noexcept {
    const double max = std::numeric_limits<double>::max();
    if (left > max - right) {
        return max;
    }
    return left + right;
}

void add_operations(std::size_t& counter, std::size_t amount = 1) noexcept {
    const std::size_t max = std::numeric_limits<std::size_t>::max();
    counter = amount > max - counter ? max : counter + amount;
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

void append_diagnostic(std::vector<SolverDiagnostic>& diagnostics,
                       SolverStatus status,
                       std::string code,
                       std::string reason,
                       std::string phase,
                       std::string message,
                       const SolverMetrics& metrics,
                       std::size_t row = SolverDiagnostic::invalid_index(),
                       std::size_t column = SolverDiagnostic::invalid_index(),
                       double value = 0.0,
                       double tolerance = 0.0) {
    diagnostics.push_back(diagnostic(status,
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

MatrixSolverResult make_matrix_result(SolverStatus status,
                                      std::string code,
                                      std::string reason,
                                      std::string phase,
                                      std::string message,
                                      SolverMetrics metrics = {},
                                      std::size_t row = SolverDiagnostic::invalid_index(),
                                      std::size_t column = SolverDiagnostic::invalid_index(),
                                      double value = 0.0,
                                      double tolerance = 0.0) {
    MatrixSolverResult result;
    result.status = status;
    result.method = "dense_partial_pivot";
    result.implementation = "MatCal::Linalg pivoted LU";
    result.metrics = metrics;
    if (!message.empty()) {
        append_diagnostic(result.diagnostics,
                          status,
                          std::move(code),
                          std::move(reason),
                          std::move(phase),
                          std::move(message),
                          result.metrics,
                          row,
                          column,
                          value,
                          tolerance);
    }
    return result;
}

PivotedLuFactorizationResult make_factorization_result(
    SolverStatus status,
    std::string code,
    std::string reason,
    std::string phase,
    std::string message,
    SolverMetrics metrics = {},
    std::size_t row = SolverDiagnostic::invalid_index(),
    std::size_t column = SolverDiagnostic::invalid_index(),
    double value = 0.0,
    double tolerance = 0.0) {
    PivotedLuFactorizationResult result;
    result.status = status;
    result.method = "dense_partial_pivot_lu";
    result.implementation = "MatCal::Linalg pivoted LU";
    result.metrics = metrics;
    if (!message.empty()) {
        append_diagnostic(result.diagnostics,
                          status,
                          std::move(code),
                          std::move(reason),
                          std::move(phase),
                          std::move(message),
                          result.metrics,
                          row,
                          column,
                          value,
                          tolerance);
    }
    return result;
}

MatrixSolverResult breakdown_matrix_result(std::string reason,
                                           std::string phase,
                                           SolverMetrics metrics,
                                           std::size_t row,
                                           std::size_t column,
                                           double value) {
    return make_matrix_result(SolverStatus::breakdown,
                              "non_finite_intermediate",
                              std::move(reason),
                              std::move(phase),
                              "Dense solve produced a non-finite intermediate value",
                              metrics,
                              row,
                              column,
                              value);
}

PivotedLuFactorizationResult breakdown_factorization_result(std::string reason,
                                                            std::string phase,
                                                            SolverMetrics metrics,
                                                            std::size_t row,
                                                            std::size_t column,
                                                            double value) {
    return make_factorization_result(SolverStatus::breakdown,
                                     "non_finite_intermediate",
                                     std::move(reason),
                                     std::move(phase),
                                     "Dense LU factorization produced a non-finite intermediate value",
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

double residual_norm_inf_checked(const DenseMatrix& a,
                                 const DenseMatrix& x,
                                 const DenseMatrix& b,
                                 SolverMetrics& metrics,
                                 MatrixSolverResult& failure) {
    double max_residual = 0.0;
    for (DenseMatrix::size_type r = 0; r < a.rows(); ++r) {
        for (DenseMatrix::size_type rhs_col = 0; rhs_col < b.cols(); ++rhs_col) {
            double sum = 0.0;
            for (DenseMatrix::size_type c = 0; c < a.cols(); ++c) {
                sum += a(r, c) * x(c, rhs_col);
                add_operations(metrics.operation_count);
                add_operations(metrics.solve_operation_count);
                if (!std::isfinite(sum)) {
                    failure = breakdown_matrix_result("residual_sum_not_finite", "residual", metrics, r, c, sum);
                    return std::numeric_limits<double>::quiet_NaN();
                }
            }
            const double residual = sum - b(r, rhs_col);
            if (!std::isfinite(residual)) {
                failure = breakdown_matrix_result("residual_not_finite",
                                                  "residual",
                                                  metrics,
                                                  r,
                                                  rhs_col,
                                                  residual);
                return std::numeric_limits<double>::quiet_NaN();
            }
            max_residual = std::max(max_residual, std::abs(residual));
        }
    }
    return max_residual;
}

} // namespace

struct PivotedLuBuilder {
    static PivotedLuFactorization create(DenseMatrix original,
                                         DenseMatrix lu,
                                         std::vector<DenseMatrix::size_type> permutation,
                                         int permutation_sign,
                                         SolverMetrics metrics) {
        PivotedLuFactorization factorization;
        factorization.original_ = std::move(original);
        factorization.lu_ = std::move(lu);
        factorization.permutation_ = std::move(permutation);
        factorization.permutation_sign_ = permutation_sign;
        factorization.metrics_ = metrics;
        factorization.valid_ = true;
        return factorization;
    }
};

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

PivotedLuFactorization::size_type PivotedLuFactorization::size() const noexcept {
    return lu_.rows();
}

bool PivotedLuFactorization::valid() const noexcept {
    return valid_;
}

const DenseMatrix& PivotedLuFactorization::lu() const noexcept {
    return lu_;
}

const DenseMatrix& PivotedLuFactorization::original_matrix() const noexcept {
    return original_;
}

const std::vector<PivotedLuFactorization::size_type>&
PivotedLuFactorization::row_permutation() const noexcept {
    return permutation_;
}

int PivotedLuFactorization::permutation_sign() const noexcept {
    return permutation_sign_;
}

const SolverMetrics& PivotedLuFactorization::factorization_metrics() const noexcept {
    return metrics_;
}

double PivotedLuFactorization::determinant() const {
    if (!valid_) {
        throw std::logic_error("Pivoted LU factorization is not valid");
    }
    double determinant_value = static_cast<double>(permutation_sign_);
    for (size_type i = 0; i < size(); ++i) {
        determinant_value *= lu_(i, i);
    }
    return determinant_value;
}

MatrixSolverResult PivotedLuFactorization::solve(const DenseMatrix& rhs,
                                                 const SolverOptions& options) const {
    if (!valid_) {
        return make_matrix_result(SolverStatus::invalid_input,
                                  "invalid_factorization",
                                  "factorization_not_valid",
                                  "validation",
                                  "Pivoted LU factorization is not valid");
    }
    if (!options.valid()) {
        return make_matrix_result(SolverStatus::invalid_input,
                                  "invalid_options",
                                  "invalid_solver_options",
                                  "validation",
                                  "Invalid solver options",
                                  metrics_);
    }
    if (rhs.rows() != size()) {
        return make_matrix_result(SolverStatus::dimension_mismatch,
                                  "rhs_size_mismatch",
                                  "dimension_mismatch",
                                  "validation",
                                  "Right-hand side row count does not match factorization size",
                                  metrics_);
    }
    if (!rhs.all_finite()) {
        return make_matrix_result(SolverStatus::non_finite_input,
                                  "non_finite_rhs",
                                  "input_not_finite",
                                  "validation",
                                  "Dense solve right-hand side must be finite",
                                  metrics_);
    }

    SolverMetrics metrics = metrics_;
    metrics.rhs_scale = max_abs_matrix_scale(rhs);
    metrics.solution_scale = 0.0;
    metrics.absolute_residual_norm = 0.0;
    metrics.relative_residual_norm = 0.0;
    metrics.residual_norm = 0.0;
    metrics.residual_acceptance_tolerance = 0.0;
    metrics.solve_operation_count = 0;
    metrics.operation_count = metrics.factorization_operation_count;

    if (!std::isfinite(metrics.rhs_scale)) {
        return breakdown_matrix_result("rhs_scale_not_finite",
                                       "scale",
                                       metrics,
                                       SolverDiagnostic::invalid_index(),
                                       SolverDiagnostic::invalid_index(),
                                       metrics.rhs_scale);
    }

    const auto n = size();
    DenseMatrix solution(n, rhs.cols());
    if (n == 0 || rhs.cols() == 0) {
        MatrixSolverResult result;
        result.status = SolverStatus::success;
        result.method = "dense_partial_pivot";
        result.implementation = "MatCal::Linalg pivoted LU";
        result.solution = std::move(solution);
        result.metrics = metrics;
        return result;
    }

    DenseMatrix y(n, rhs.cols());
    for (size_type r = 0; r < n; ++r) {
        for (size_type col = 0; col < rhs.cols(); ++col) {
            y(r, col) = rhs(permutation_[r], col);
        }
    }

    for (size_type i = 0; i < n; ++i) {
        for (size_type col = 0; col < rhs.cols(); ++col) {
            double sum = y(i, col);
            for (size_type k = 0; k < i; ++k) {
                sum -= lu_(i, k) * y(k, col);
                add_operations(metrics.operation_count);
                add_operations(metrics.solve_operation_count);
                if (!std::isfinite(sum)) {
                    return breakdown_matrix_result("forward_substitution_sum_not_finite",
                                                   "forward_substitution",
                                                   metrics,
                                                   i,
                                                   k,
                                                   sum);
                }
            }
            y(i, col) = sum;
        }
    }

    for (size_type i = n; i-- > 0;) {
        const double diag = lu_(i, i);
        const double diag_abs = std::abs(diag);
        if (!std::isfinite(diag_abs)) {
            return breakdown_matrix_result("back_substitution_pivot_not_finite",
                                           "back_substitution",
                                           metrics,
                                           i,
                                           i,
                                           diag);
        }
        if (diag_abs <= metrics.pivot_tolerance_used) {
            return make_matrix_result(SolverStatus::singular,
                                      "pivot_too_small",
                                      "back_substitution_pivot_below_tolerance",
                                      "back_substitution",
                                      "Back substitution pivot is below tolerance",
                                      metrics,
                                      i,
                                      i,
                                      diag_abs,
                                      metrics.pivot_tolerance_used);
        }
        for (size_type col = 0; col < rhs.cols(); ++col) {
            double sum = y(i, col);
            for (size_type k = i + 1; k < n; ++k) {
                sum -= lu_(i, k) * solution(k, col);
                add_operations(metrics.operation_count);
                add_operations(metrics.solve_operation_count);
                if (!std::isfinite(sum)) {
                    return breakdown_matrix_result("back_substitution_sum_not_finite",
                                                   "back_substitution",
                                                   metrics,
                                                   i,
                                                   k,
                                                   sum);
                }
            }
            solution(i, col) = sum / diag;
            if (!std::isfinite(solution(i, col))) {
                return breakdown_matrix_result("solution_not_finite",
                                               "back_substitution",
                                               metrics,
                                               i,
                                               col,
                                               solution(i, col));
            }
        }
    }

    metrics.solution_scale = max_abs_matrix_scale(solution);
    if (!std::isfinite(metrics.solution_scale)) {
        return breakdown_matrix_result("solution_scale_not_finite",
                                       "solution",
                                       metrics,
                                       SolverDiagnostic::invalid_index(),
                                       SolverDiagnostic::invalid_index(),
                                       metrics.solution_scale);
    }

    MatrixSolverResult residual_failure;
    const double max_residual = residual_norm_inf_checked(original_, solution, rhs, metrics, residual_failure);
    if (!std::isfinite(max_residual)) {
        return residual_failure;
    }

    metrics.absolute_residual_norm = max_residual;
    metrics.residual_norm = max_residual;

    const double denominator = saturated_sum(saturated_product(metrics.matrix_scale, metrics.solution_scale),
                                             metrics.rhs_scale);
    if (denominator == 0.0) {
        metrics.relative_residual_norm = max_residual == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    } else {
        metrics.relative_residual_norm = max_residual / denominator;
    }
    if (!std::isfinite(metrics.relative_residual_norm)) {
        return breakdown_matrix_result("relative_residual_not_finite",
                                       "residual",
                                       metrics,
                                       SolverDiagnostic::invalid_index(),
                                       SolverDiagnostic::invalid_index(),
                                       metrics.relative_residual_norm);
    }

    metrics.residual_acceptance_tolerance = options.comparison_tolerance(denominator);
    if (!std::isfinite(metrics.residual_acceptance_tolerance)) {
        return breakdown_matrix_result("residual_tolerance_not_finite",
                                       "residual",
                                       metrics,
                                       SolverDiagnostic::invalid_index(),
                                       SolverDiagnostic::invalid_index(),
                                       metrics.residual_acceptance_tolerance);
    }
    if (metrics.absolute_residual_norm > metrics.residual_acceptance_tolerance) {
        return make_matrix_result(SolverStatus::not_converged,
                                  "residual_too_large",
                                  "residual_above_tolerance",
                                  "residual",
                                  "Dense solve residual is above acceptance tolerance",
                                  metrics,
                                  SolverDiagnostic::invalid_index(),
                                  SolverDiagnostic::invalid_index(),
                                  metrics.absolute_residual_norm,
                                  metrics.residual_acceptance_tolerance);
    }

    MatrixSolverResult result;
    result.status = SolverStatus::success;
    result.method = "dense_partial_pivot";
    result.implementation = "MatCal::Linalg pivoted LU";
    result.solution = std::move(solution);
    result.metrics = metrics;
    return result;
}

SolverResult PivotedLuFactorization::solve(const Vector& rhs,
                                           const SolverOptions& options) const {
    DenseMatrix rhs_matrix(rhs.size(), 1);
    for (size_type i = 0; i < rhs.size(); ++i) {
        rhs_matrix(i, 0) = rhs[i];
    }

    auto matrix_result = solve(rhs_matrix, options);
    SolverResult result;
    result.status = matrix_result.status;
    result.metrics = matrix_result.metrics;
    result.diagnostics = std::move(matrix_result.diagnostics);
    result.method = std::move(matrix_result.method);
    result.implementation = std::move(matrix_result.implementation);
    if (matrix_result.success()) {
        result.solution = Vector(rhs.size());
        for (size_type i = 0; i < rhs.size(); ++i) {
            result.solution[i] = matrix_result.solution(i, 0);
        }
    }
    return result;
}

PivotedLuFactorizationResult factorize_dense_partial_pivot(const DenseMatrix& a,
                                                           const SolverOptions& options) {
    if (!options.valid()) {
        return make_factorization_result(SolverStatus::invalid_input,
                                         "invalid_options",
                                         "invalid_solver_options",
                                         "validation",
                                         "Invalid solver options");
    }
    if (a.rows() != a.cols()) {
        return make_factorization_result(SolverStatus::dimension_mismatch,
                                         "matrix_not_square",
                                         "dimension_mismatch",
                                         "validation",
                                         "Dense LU factorization requires a square matrix");
    }
    if (!a.all_finite()) {
        return make_factorization_result(SolverStatus::non_finite_input,
                                         "non_finite_matrix",
                                         "input_not_finite",
                                         "validation",
                                         "Dense LU factorization input must be finite");
    }

    const auto n = a.rows();
    SolverMetrics metrics;
    metrics.matrix_scale = max_abs_matrix_scale(a);
    metrics.pivot_tolerance_used = options.pivot_tolerance(metrics.matrix_scale);
    metrics.minimum_abs_pivot = n == 0 ? 0.0 : std::numeric_limits<double>::infinity();

    if (!std::isfinite(metrics.matrix_scale) || !std::isfinite(metrics.pivot_tolerance_used)) {
        return breakdown_factorization_result("scale_not_finite",
                                              "scale",
                                              metrics,
                                              SolverDiagnostic::invalid_index(),
                                              SolverDiagnostic::invalid_index(),
                                              metrics.matrix_scale);
    }

    std::vector<double> work = a.values();
    std::vector<DenseMatrix::size_type> permutation(n);
    for (DenseMatrix::size_type i = 0; i < n; ++i) {
        permutation[i] = i;
    }
    int permutation_sign = 1;

    for (DenseMatrix::size_type k = 0; k < n; ++k) {
        DenseMatrix::size_type pivot_row = k;
        double pivot_abs = std::abs(work[k * n + k]);
        for (DenseMatrix::size_type r = k + 1; r < n; ++r) {
            const double candidate = std::abs(work[r * n + k]);
            if (!std::isfinite(candidate)) {
                return breakdown_factorization_result("pivot_search_not_finite",
                                                      "pivot_search",
                                                      metrics,
                                                      r,
                                                      k,
                                                      candidate);
            }
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot_row = r;
            }
        }

        if (!std::isfinite(pivot_abs)) {
            return breakdown_factorization_result("pivot_not_finite",
                                                  "pivot_search",
                                                  metrics,
                                                  pivot_row,
                                                  k,
                                                  pivot_abs);
        }
        metrics.minimum_abs_pivot = std::min(metrics.minimum_abs_pivot, pivot_abs);

        if (pivot_abs <= metrics.pivot_tolerance_used) {
            return make_factorization_result(SolverStatus::singular,
                                             "pivot_too_small",
                                             "pivot_below_tolerance",
                                             "factorization",
                                             "Pivot is zero or below pivot tolerance",
                                             metrics,
                                             pivot_row,
                                             k,
                                             pivot_abs,
                                             metrics.pivot_tolerance_used);
        }

        if (pivot_row != k) {
            swap_rows(work, n, k, pivot_row);
            std::swap(permutation[k], permutation[pivot_row]);
            permutation_sign = -permutation_sign;
        }

        const double pivot = work[k * n + k];
        if (!std::isfinite(pivot)) {
            return breakdown_factorization_result("pivot_not_finite", "factorization", metrics, k, k, pivot);
        }
        for (DenseMatrix::size_type r = k + 1; r < n; ++r) {
            const double factor = work[r * n + k] / pivot;
            if (!std::isfinite(factor)) {
                return breakdown_factorization_result("factor_not_finite", "elimination", metrics, r, k, factor);
            }
            work[r * n + k] = factor;
            for (DenseMatrix::size_type c = k + 1; c < n; ++c) {
                work[r * n + c] -= factor * work[k * n + c];
                add_operations(metrics.operation_count);
                if (!std::isfinite(work[r * n + c])) {
                    return breakdown_factorization_result("matrix_update_not_finite",
                                                          "elimination",
                                                          metrics,
                                                          r,
                                                          c,
                                                          work[r * n + c]);
                }
            }
        }
        metrics.iterations = k + 1;
    }
    metrics.factorization_operation_count = metrics.operation_count;

    auto lu = DenseMatrix::from_row_major(n, n, std::move(work));
    PivotedLuFactorizationResult result;
    result.status = SolverStatus::success;
    result.factorization = PivotedLuBuilder::create(a, std::move(lu), std::move(permutation), permutation_sign, metrics);
    result.metrics = metrics;
    result.method = "dense_partial_pivot_lu";
    result.implementation = "MatCal::Linalg pivoted LU";
    return result;
}

SolverResult solve_dense_partial_pivot(const DenseMatrix& a,
                                       const Vector& b,
                                       const SolverOptions& options) {
    auto factorization = factorize_dense_partial_pivot(a, options);
    if (!factorization.success()) {
        SolverResult result;
        result.status = factorization.status;
        result.metrics = factorization.metrics;
        result.diagnostics = std::move(factorization.diagnostics);
        result.method = "dense_partial_pivot";
        result.implementation = "MatCal::Linalg pivoted LU";
        return result;
    }
    return factorization.factorization.solve(b, options);
}

} // namespace MatCal::Linalg
