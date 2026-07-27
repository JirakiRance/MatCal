#include "MatCal/Linalg/SkylineLdlt.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

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

double vector_scale(const Vector& vector) noexcept {
    double scale = 0.0;
    for (double value : vector) {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

bool add_operations(SolverMetrics& metrics, std::size_t amount, bool solve_phase) noexcept {
    if (amount > std::numeric_limits<std::size_t>::max() - metrics.operation_count) {
        return false;
    }
    if (solve_phase) {
        if (amount > std::numeric_limits<std::size_t>::max() - metrics.solve_operation_count) {
            return false;
        }
    } else {
        if (amount > std::numeric_limits<std::size_t>::max() - metrics.factorization_operation_count) {
            return false;
        }
    }

    metrics.operation_count += amount;
    if (solve_phase) {
        metrics.solve_operation_count += amount;
    } else {
        metrics.factorization_operation_count += amount;
    }
    return true;
}

SolverDiagnostic make_diagnostic(SolverStatus status,
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

SkylineLdltFactorizationResult factor_result(SolverStatus status,
                                             std::string code,
                                             std::string reason,
                                             std::string phase,
                                             std::string message,
                                             SolverMetrics metrics = {},
                                             std::size_t row = SolverDiagnostic::invalid_index(),
                                             std::size_t column = SolverDiagnostic::invalid_index(),
                                             double value = 0.0,
                                             double tolerance = 0.0) {
    SkylineLdltFactorizationResult result;
    result.status = status;
    result.metrics = metrics;
    if (!message.empty()) {
        result.diagnostics.push_back(make_diagnostic(status,
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

SolverResult solve_result(SolverStatus status,
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
    result.method = "skyline_ldlt";
    result.implementation = "MatCal::Linalg skyline reference";
    result.metrics = metrics;
    if (!message.empty()) {
        result.diagnostics.push_back(make_diagnostic(status,
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

SolverResult breakdown_solve(std::string reason,
                             std::string phase,
                             SolverMetrics metrics,
                             std::size_t row,
                             std::size_t column,
                             double value) {
    return solve_result(SolverStatus::breakdown,
                        "non_finite_intermediate",
                        std::move(reason),
                        std::move(phase),
                        "Skyline LDLT solve produced a non-finite intermediate value",
                        metrics,
                        row,
                        column,
                        value);
}

SolverResult operation_overflow_solve(std::string phase, SolverMetrics metrics) {
    return solve_result(SolverStatus::breakdown,
                        "operation_count_overflow",
                        "operation_count_overflow",
                        std::move(phase),
                        "Skyline LDLT solve operation count overflowed",
                        metrics);
}

SkylineLdltFactorizationResult operation_overflow_factor(SolverMetrics metrics) {
    return factor_result(SolverStatus::breakdown,
                         "operation_count_overflow",
                         "operation_count_overflow",
                         "factorization",
                         "Skyline LDLT factorization operation count overflowed",
                         metrics);
}

} // namespace

double SkylineLdltFactorization::diagonal(std::size_t index) const {
    if (index >= size_) {
        throw std::out_of_range("LDLT diagonal index out of range");
    }
    return diagonal_[index];
}

double SkylineLdltFactorization::lower(std::size_t row, std::size_t column) const {
    if (row >= size_ || column >= size_) {
        throw std::out_of_range("LDLT lower index out of range");
    }
    if (row == column) {
        return 1.0;
    }
    if (row < column || !stores(row, column)) {
        return 0.0;
    }
    return l_values_[storage_index(row, column)];
}

Vector SkylineLdltFactorization::multiply(const Vector& x) const {
    if (x.size() != size_) {
        throw std::invalid_argument("LDLT matvec dimension mismatch");
    }

    Vector transposed(size_);
    for (std::size_t i = 0; i < size_; ++i) {
        transposed[i] = x[i];
    }
    for (std::size_t row = 0; row < size_; ++row) {
        for (std::size_t column = first_columns_[row]; column < row; ++column) {
            transposed[column] += lower(row, column) * x[row];
        }
    }

    for (std::size_t i = 0; i < size_; ++i) {
        transposed[i] *= diagonal_[i];
    }

    Vector result(size_);
    for (std::size_t row = 0; row < size_; ++row) {
        double sum = transposed[row];
        for (std::size_t column = first_columns_[row]; column < row; ++column) {
            sum += lower(row, column) * transposed[column];
        }
        result[row] = sum;
    }
    return result;
}

SolverResult SkylineLdltFactorization::solve(const Vector& b, const SolverOptions& options) const {
    SolverMetrics metrics;
    metrics.matrix_scale = matrix_scale_;
    metrics.rhs_scale = vector_scale(b);
    metrics.pivot_tolerance_used = pivot_tolerance_used_;
    metrics.minimum_abs_pivot = minimum_abs_pivot_;
    metrics.factorization_operation_count = factorization_operation_count_;
    metrics.operation_count = factorization_operation_count_;
    metrics.iterations = size_;

    if (!options.valid()) {
        return solve_result(SolverStatus::invalid_input,
                            "invalid_options",
                            "invalid_solver_options",
                            "validation",
                            "Invalid solver options",
                            metrics);
    }
    if (b.size() != size_) {
        return solve_result(SolverStatus::dimension_mismatch,
                            "rhs_size_mismatch",
                            "dimension_mismatch",
                            "validation",
                            "Right-hand side size does not match factorization size",
                            metrics);
    }
    if (!b.all_finite()) {
        return solve_result(SolverStatus::non_finite_input,
                            "non_finite_rhs",
                            "rhs_not_finite",
                            "validation",
                            "Skyline LDLT right-hand side must be finite",
                            metrics);
    }
    if (!std::isfinite(metrics.rhs_scale)) {
        return breakdown_solve("rhs_scale_not_finite", "scale", metrics, SolverDiagnostic::invalid_index(), SolverDiagnostic::invalid_index(), metrics.rhs_scale);
    }

    Vector y(size_);
    for (std::size_t i = 0; i < size_; ++i) {
        double sum = b[i];
        for (std::size_t j = first_columns_[i]; j < i; ++j) {
            sum -= lower(i, j) * y[j];
            if (!add_operations(metrics, 1, true)) {
                return operation_overflow_solve("forward_substitution", metrics);
            }
            if (!std::isfinite(sum)) {
                return breakdown_solve("forward_substitution_not_finite", "forward_substitution", metrics, i, j, sum);
            }
        }
        y[i] = sum;
    }

    Vector z(size_);
    for (std::size_t i = 0; i < size_; ++i) {
        z[i] = y[i] / diagonal_[i];
        if (!add_operations(metrics, 1, true)) {
            return operation_overflow_solve("diagonal_solve", metrics);
        }
        if (!std::isfinite(z[i])) {
            return breakdown_solve("diagonal_solve_not_finite", "diagonal_solve", metrics, i, i, z[i]);
        }
    }

    Vector x(size_);
    for (std::size_t i = size_; i-- > 0;) {
        double sum = z[i];
        for (std::size_t row : lower_rows_by_column_[i]) {
            sum -= lower(row, i) * x[row];
            if (!add_operations(metrics, 1, true)) {
                return operation_overflow_solve("back_substitution", metrics);
            }
            if (!std::isfinite(sum)) {
                return breakdown_solve("back_substitution_not_finite", "back_substitution", metrics, row, i, sum);
            }
        }
        x[i] = sum;
    }

    metrics.solution_scale = vector_scale(x);
    if (!std::isfinite(metrics.solution_scale)) {
        return breakdown_solve("solution_scale_not_finite", "solution", metrics, SolverDiagnostic::invalid_index(), SolverDiagnostic::invalid_index(), metrics.solution_scale);
    }

    Vector ax = multiply(x);
    if (!add_operations(metrics, storage_size(), true)) {
        return operation_overflow_solve("residual", metrics);
    }
    double absolute_residual = 0.0;
    for (std::size_t i = 0; i < size_; ++i) {
        double residual = ax[i] - b[i];
        if (!std::isfinite(residual)) {
            return breakdown_solve("residual_not_finite", "residual", metrics, i, SolverDiagnostic::invalid_index(), residual);
        }
        absolute_residual = std::max(absolute_residual, std::abs(residual));
    }

    metrics.absolute_residual_norm = absolute_residual;
    metrics.residual_norm = absolute_residual;
    double denominator = saturated_sum(saturated_product(metrics.matrix_scale, metrics.solution_scale),
                                      metrics.rhs_scale);
    metrics.relative_residual_norm = denominator == 0.0 ? (absolute_residual == 0.0 ? 0.0 : std::numeric_limits<double>::infinity())
                                                        : absolute_residual / denominator;
    if (!std::isfinite(metrics.relative_residual_norm)) {
        return breakdown_solve("relative_residual_not_finite", "residual", metrics, SolverDiagnostic::invalid_index(), SolverDiagnostic::invalid_index(), metrics.relative_residual_norm);
    }

    metrics.residual_acceptance_tolerance = options.comparison_tolerance(denominator);
    if (!std::isfinite(metrics.residual_acceptance_tolerance)) {
        return breakdown_solve("residual_tolerance_not_finite", "residual", metrics, SolverDiagnostic::invalid_index(), SolverDiagnostic::invalid_index(), metrics.residual_acceptance_tolerance);
    }
    if (metrics.absolute_residual_norm > metrics.residual_acceptance_tolerance) {
        return solve_result(SolverStatus::not_converged,
                            "residual_too_large",
                            "residual_above_tolerance",
                            "residual",
                            "Skyline LDLT residual is above acceptance tolerance",
                            metrics,
                            SolverDiagnostic::invalid_index(),
                            SolverDiagnostic::invalid_index(),
                            metrics.absolute_residual_norm,
                            metrics.residual_acceptance_tolerance);
    }

    SolverResult result;
    result.status = SolverStatus::success;
    result.solution = std::move(x);
    result.metrics = metrics;
    result.method = "skyline_ldlt";
    result.implementation = "MatCal::Linalg skyline reference";
    return result;
}

bool SkylineLdltFactorization::stores(std::size_t row, std::size_t column) const {
    return row < size_ && column < size_ && row >= column && column >= first_columns_[row];
}

std::size_t SkylineLdltFactorization::storage_index(std::size_t row, std::size_t column) const {
    return row_offsets_[row] + (column - first_columns_[row]);
}

SkylineLdltFactorizationResult factorize_skyline_ldlt(const SymmetricSkylineMatrix& matrix,
                                                       const SolverOptions& options) {
    SolverMetrics metrics;
    metrics.matrix_scale = matrix.matrix_scale();
    metrics.pivot_tolerance_used = options.pivot_tolerance(metrics.matrix_scale);

    if (!options.valid()) {
        return factor_result(SolverStatus::invalid_input,
                             "invalid_options",
                             "invalid_solver_options",
                             "validation",
                             "Invalid solver options",
                             metrics);
    }
    if (!matrix.all_finite()) {
        return factor_result(SolverStatus::non_finite_input,
                             "non_finite_matrix",
                             "matrix_not_finite",
                             "validation",
                             "Skyline LDLT input matrix must be finite",
                             metrics);
    }
    if (!std::isfinite(metrics.matrix_scale) || !std::isfinite(metrics.pivot_tolerance_used)) {
        return factor_result(SolverStatus::breakdown,
                             "non_finite_intermediate",
                             "scale_not_finite",
                             "scale",
                             "Skyline LDLT scale computation produced a non-finite value",
                             metrics,
                             SolverDiagnostic::invalid_index(),
                             SolverDiagnostic::invalid_index(),
                             metrics.matrix_scale);
    }

    SkylineLdltFactorization factorization;
    factorization.size_ = matrix.size();
    factorization.first_columns_ = matrix.first_columns();
    factorization.row_offsets_ = matrix.row_offsets();
    factorization.lower_rows_by_column_.assign(matrix.size(), {});
    factorization.l_values_.assign(matrix.storage_size(), 0.0);
    factorization.diagonal_.assign(matrix.size(), 0.0);
    factorization.matrix_scale_ = metrics.matrix_scale;
    factorization.pivot_tolerance_used_ = metrics.pivot_tolerance_used;
    factorization.minimum_abs_pivot_ = matrix.size() == 0 ? 0.0 : std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < matrix.size(); ++i) {
        factorization.l_values_[factorization.storage_index(i, i)] = 1.0;
        for (std::size_t column = matrix.first_columns()[i]; column < i; ++column) {
            factorization.lower_rows_by_column_[column].push_back(i);
        }
    }

    for (std::size_t i = 0; i < matrix.size(); ++i) {
        for (std::size_t j = matrix.first_columns()[i]; j < i; ++j) {
            double sum = matrix.get(i, j);
            std::size_t k_start = std::max(matrix.first_columns()[i], matrix.first_columns()[j]);
            for (std::size_t k = k_start; k < j; ++k) {
                double term = factorization.lower(i, k) * factorization.diagonal_[k] * factorization.lower(j, k);
                sum -= term;
                if (!add_operations(metrics, 1, false)) {
                    return operation_overflow_factor(metrics);
                }
                if (!std::isfinite(sum)) {
                    return factor_result(SolverStatus::breakdown,
                                         "non_finite_intermediate",
                                         "offdiagonal_update_not_finite",
                                         "factorization",
                                         "Skyline LDLT offdiagonal update produced a non-finite value",
                                         metrics,
                                         i,
                                         j,
                                         sum);
                }
            }
            double pivot = factorization.diagonal_[j];
            double value = sum / pivot;
            if (!std::isfinite(value)) {
                return factor_result(SolverStatus::breakdown,
                                     "non_finite_intermediate",
                                     "offdiagonal_factor_not_finite",
                                     "factorization",
                                     "Skyline LDLT offdiagonal factor produced a non-finite value",
                                     metrics,
                                     i,
                                     j,
                                     value);
            }
            factorization.l_values_[factorization.storage_index(i, j)] = value;
        }

        double pivot = matrix.get(i, i);
        for (std::size_t k = matrix.first_columns()[i]; k < i; ++k) {
            double l = factorization.lower(i, k);
            pivot -= l * l * factorization.diagonal_[k];
            if (!add_operations(metrics, 1, false)) {
                return operation_overflow_factor(metrics);
            }
            if (!std::isfinite(pivot)) {
                return factor_result(SolverStatus::breakdown,
                                     "non_finite_intermediate",
                                     "diagonal_update_not_finite",
                                     "factorization",
                                     "Skyline LDLT diagonal update produced a non-finite value",
                                     metrics,
                                     i,
                                     k,
                                     pivot);
            }
        }

        double abs_pivot = std::abs(pivot);
        metrics.minimum_abs_pivot = std::min(metrics.minimum_abs_pivot, abs_pivot);
        factorization.minimum_abs_pivot_ = std::min(factorization.minimum_abs_pivot_, abs_pivot);
        if (!std::isfinite(abs_pivot)) {
            return factor_result(SolverStatus::breakdown,
                                 "non_finite_intermediate",
                                 "pivot_not_finite",
                                 "factorization",
                                 "Skyline LDLT pivot is not finite",
                                 metrics,
                                 i,
                                 i,
                                 pivot);
        }
        if (pivot <= metrics.pivot_tolerance_used) {
            return factor_result(SolverStatus::not_positive_definite,
                                 "non_positive_pivot",
                                 "pivot_not_positive_definite",
                                 "factorization",
                                 "Skyline LDLT pivot is not positive above tolerance",
                                 metrics,
                                 i,
                                 i,
                                 pivot,
                                 metrics.pivot_tolerance_used);
        }

        factorization.diagonal_[i] = pivot;
        metrics.iterations = i + 1;
    }

    if (matrix.size() == 0) {
        metrics.minimum_abs_pivot = 0.0;
    }
    factorization.factorization_operation_count_ = metrics.operation_count;
    metrics.factorization_operation_count = metrics.operation_count;

    SkylineLdltFactorizationResult result;
    result.status = SolverStatus::success;
    result.factorization = std::move(factorization);
    result.metrics = metrics;
    return result;
}

SolverResult solve_skyline_ldlt(const SymmetricSkylineMatrix& matrix,
                                const Vector& b,
                                const SolverOptions& options) {
    auto factorized = factorize_skyline_ldlt(matrix, options);
    if (!factorized.success()) {
        SolverResult result;
        result.status = factorized.status;
        result.metrics = factorized.metrics;
        result.diagnostics = factorized.diagnostics;
        result.method = "skyline_ldlt";
        result.implementation = "MatCal::Linalg skyline reference";
        return result;
    }
    return factorized.factorization.solve(b, options);
}

} // namespace MatCal::Linalg
