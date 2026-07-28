#ifndef MATCAL_LINALG_EIGEN_SOLVERS_HPP
#define MATCAL_LINALG_EIGEN_SOLVERS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/SolverTypes.hpp"
#include "MatCal/Linalg/Vector.hpp"

namespace MatCal::Linalg {

enum class EigenStatus {
    success,
    invalid_input,
    dimension_mismatch,
    non_finite_input,
    singular_shift,
    breakdown,
    not_converged
};

inline const char* to_string(EigenStatus status) noexcept {
    switch (status) {
    case EigenStatus::success: return "success";
    case EigenStatus::invalid_input: return "invalid_input";
    case EigenStatus::dimension_mismatch: return "dimension_mismatch";
    case EigenStatus::non_finite_input: return "non_finite_input";
    case EigenStatus::singular_shift: return "singular_shift";
    case EigenStatus::breakdown: return "breakdown";
    case EigenStatus::not_converged: return "not_converged";
    }
    return "unknown";
}

struct EigenOptions {
    double absolute_tolerance = 0.0;
    double relative_tolerance = 1.0e-10;
    std::size_t max_iterations = 1000;
    double shift = 0.0;
    SolverOptions linear_options{};

    bool valid() const noexcept {
        return std::isfinite(absolute_tolerance) &&
               std::isfinite(relative_tolerance) &&
               std::isfinite(shift) &&
               absolute_tolerance >= 0.0 &&
               relative_tolerance >= 0.0 &&
               max_iterations > 0 &&
               linear_options.valid();
    }

    double comparison_tolerance(double scale) const noexcept {
        if (!std::isfinite(scale) || scale < 0.0) {
            return std::numeric_limits<double>::max();
        }
        if (relative_tolerance != 0.0 && scale > std::numeric_limits<double>::max() / relative_tolerance) {
            return std::numeric_limits<double>::max();
        }
        return std::max(absolute_tolerance, relative_tolerance * scale);
    }
};

struct EigenDiagnostic {
    EigenStatus status = EigenStatus::success;
    std::string code;
    std::string reason;
    std::string phase;
    std::size_t row = invalid_index();
    std::size_t column = invalid_index();
    double value = 0.0;
    std::string message;

    static constexpr std::size_t invalid_index() noexcept {
        return std::numeric_limits<std::size_t>::max();
    }
};

struct EigenMetrics {
    std::size_t iterations = 0;
    std::size_t operation_count = 0;
    std::size_t linear_solves = 0;
    double residual_norm = std::numeric_limits<double>::infinity();
    double residual_acceptance_tolerance = 0.0;
    double matrix_scale = 0.0;
    double eigenvector_scale = 0.0;
    double eigenvalue_delta = std::numeric_limits<double>::infinity();
    double shift = 0.0;
};

struct EigenResult {
    EigenStatus status = EigenStatus::invalid_input;
    double eigenvalue = 0.0;
    Vector eigenvector;
    EigenMetrics metrics;
    std::vector<EigenDiagnostic> diagnostics;
    std::string method;

    bool success() const noexcept {
        return status == EigenStatus::success;
    }
};

inline EigenResult dominant_eigenpair(const DenseMatrix& a,
                                      const EigenOptions& options = {});
inline EigenResult dominant_eigenpair(const DenseMatrix& a,
                                      const Vector& initial_vector,
                                      const EigenOptions& options = {});
inline EigenResult inverse_power_eigenpair(const DenseMatrix& a,
                                           const EigenOptions& options = {});
inline EigenResult inverse_power_eigenpair(const DenseMatrix& a,
                                           const Vector& initial_vector,
                                           const EigenOptions& options = {});

namespace eigen_detail {

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

inline EigenDiagnostic diagnostic(EigenStatus status,
                                  std::string code,
                                  std::string reason,
                                  std::string phase,
                                  std::size_t row,
                                  std::size_t column,
                                  double value,
                                  std::string message) {
    EigenDiagnostic diag;
    diag.status = status;
    diag.code = std::move(code);
    diag.reason = std::move(reason);
    diag.phase = std::move(phase);
    diag.row = row;
    diag.column = column;
    diag.value = value;
    diag.message = std::move(message);
    return diag;
}

inline EigenResult make_result(EigenStatus status,
                               std::string method,
                               std::string code,
                               std::string reason,
                               std::string phase,
                               std::string message,
                               EigenMetrics metrics = {},
                               std::size_t row = EigenDiagnostic::invalid_index(),
                               std::size_t column = EigenDiagnostic::invalid_index(),
                               double value = 0.0) {
    EigenResult result;
    result.status = status;
    result.method = std::move(method);
    result.metrics = metrics;
    if (!message.empty()) {
        result.diagnostics.push_back(diagnostic(status, std::move(code), std::move(reason), std::move(phase),
                                                row, column, value, std::move(message)));
    }
    return result;
}

inline Vector default_initial_vector(std::size_t n) {
    return Vector(n, 1.0);
}

inline bool zero_vector(const Vector& vector) noexcept {
    return vector_scale(vector) == 0.0;
}

inline std::size_t max_abs_index(const Vector& vector) {
    std::size_t index = 0;
    double best = std::abs(vector[0]);
    for (std::size_t i = 1; i < vector.size(); ++i) {
        const double candidate = std::abs(vector[i]);
        if (candidate > best) {
            best = candidate;
            index = i;
        }
    }
    return index;
}

inline double rayleigh_quotient(const DenseMatrix& a, const Vector& x, std::size_t& operation_count) {
    Vector ax = a.multiply(x);
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        numerator += x[i] * ax[i];
        denominator += x[i] * x[i];
        operation_count += 2;
    }
    if (denominator == 0.0 || !std::isfinite(numerator) || !std::isfinite(denominator)) {
        throw std::runtime_error("Rayleigh quotient is not finite");
    }
    return numerator / denominator;
}

inline double eigen_residual(const DenseMatrix& a,
                             const Vector& x,
                             double eigenvalue,
                             std::size_t& operation_count) {
    Vector ax = a.multiply(x);
    double residual = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double value = ax[i] - eigenvalue * x[i];
        ++operation_count;
        if (!std::isfinite(value)) {
            return std::numeric_limits<double>::infinity();
        }
        residual = std::max(residual, std::abs(value));
    }
    return residual;
}

inline EigenResult validate_setup(const DenseMatrix& a,
                                  const Vector* initial,
                                  const EigenOptions& options,
                                  const char* method,
                                  EigenMetrics& metrics) {
    metrics.matrix_scale = matrix_scale(a);
    metrics.shift = options.shift;
    if (!options.valid()) {
        return make_result(EigenStatus::invalid_input, method, "invalid_options",
                           "invalid_eigen_options", "validation", "Invalid eigen solver options", metrics);
    }
    if (a.rows() != a.cols()) {
        return make_result(EigenStatus::dimension_mismatch, method, "matrix_not_square",
                           "dimension_mismatch", "validation", "Eigen solve requires a square matrix", metrics);
    }
    if (initial != nullptr && initial->size() != a.rows()) {
        return make_result(EigenStatus::dimension_mismatch, method, "initial_size_mismatch",
                           "dimension_mismatch", "validation", "Initial eigenvector size mismatch", metrics);
    }
    if (!a.all_finite() || (initial != nullptr && !initial->all_finite())) {
        return make_result(EigenStatus::non_finite_input, method, "non_finite_input",
                           "input_not_finite", "validation", "Eigen inputs must be finite", metrics);
    }
    if (a.rows() == 0) {
        return make_result(EigenStatus::invalid_input, method, "empty_matrix",
                           "invalid_input", "validation", "Eigen solve requires a non-empty matrix", metrics);
    }
    EigenResult result;
    result.status = EigenStatus::success;
    result.method = method;
    result.metrics = metrics;
    return result;
}

inline EigenResult finish_success(std::string method,
                                  double eigenvalue,
                                  Vector eigenvector,
                                  EigenMetrics metrics) {
    EigenResult result = make_result(EigenStatus::success, std::move(method), "", "", "", "", metrics);
    result.eigenvalue = eigenvalue;
    result.eigenvector = std::move(eigenvector);
    return result;
}

} // namespace eigen_detail

inline EigenResult dominant_eigenpair(const DenseMatrix& a,
                                      const EigenOptions& options) {
    return dominant_eigenpair(a, eigen_detail::default_initial_vector(a.rows()), options);
}

inline EigenResult dominant_eigenpair(const DenseMatrix& a,
                                      const Vector& initial_vector,
                                      const EigenOptions& options) {
    EigenMetrics metrics;
    auto setup = eigen_detail::validate_setup(a, &initial_vector, options, "power", metrics);
    if (!setup.success()) {
        return setup;
    }
    if (eigen_detail::zero_vector(initial_vector)) {
        return eigen_detail::make_result(EigenStatus::invalid_input, "power", "zero_initial_vector",
                                         "invalid_input", "validation", "Initial eigenvector must be nonzero", metrics);
    }

    Vector x = initial_vector;
    const double initial_scale = eigen_detail::vector_scale(x);
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] /= initial_scale;
    }

    double previous_eigenvalue = 0.0;
    double eigenvalue = 0.0;
    for (std::size_t iteration = 1; iteration <= options.max_iterations; ++iteration) {
        Vector y = a.multiply(x);
        metrics.operation_count += a.rows() * a.cols();
        if (!y.all_finite()) {
            return eigen_detail::make_result(EigenStatus::breakdown, "power", "matvec_not_finite",
                                             "non_finite_intermediate", "iteration",
                                             "Power iteration produced a non-finite matrix-vector product", metrics);
        }
        const double y_scale = eigen_detail::vector_scale(y);
        if (y_scale == 0.0 || !std::isfinite(y_scale)) {
            return eigen_detail::make_result(EigenStatus::breakdown, "power", "zero_power_vector",
                                             "zero_or_non_finite_vector", "normalization",
                                             "Power iteration produced a zero or non-finite vector", metrics);
        }
        const std::size_t pivot = eigen_detail::max_abs_index(y);
        const double signed_scale = y[pivot];
        for (std::size_t i = 0; i < y.size(); ++i) {
            y[i] /= signed_scale;
        }
        x = std::move(y);

        try {
            eigenvalue = eigen_detail::rayleigh_quotient(a, x, metrics.operation_count);
            metrics.residual_norm = eigen_detail::eigen_residual(a, x, eigenvalue, metrics.operation_count);
        } catch (const std::exception& e) {
            return eigen_detail::make_result(EigenStatus::breakdown, "power", "rayleigh_breakdown",
                                             "non_finite_intermediate", "rayleigh", e.what(), metrics);
        }
        metrics.iterations = iteration;
        metrics.eigenvalue_delta = std::abs(eigenvalue - previous_eigenvalue);
        metrics.eigenvector_scale = eigen_detail::vector_scale(x);
        metrics.residual_acceptance_tolerance =
            options.comparison_tolerance(metrics.matrix_scale * metrics.eigenvector_scale +
                                         std::abs(eigenvalue) * metrics.eigenvector_scale);
        if (metrics.residual_norm <= metrics.residual_acceptance_tolerance) {
            return eigen_detail::finish_success("power", eigenvalue, x, metrics);
        }
        previous_eigenvalue = eigenvalue;
    }

    return eigen_detail::make_result(EigenStatus::not_converged, "power", "maximum_iterations",
                                     "maximum_iterations", "convergence",
                                     "Power iteration reached maximum iterations", metrics);
}

inline EigenResult inverse_power_eigenpair(const DenseMatrix& a,
                                           const EigenOptions& options) {
    return inverse_power_eigenpair(a, eigen_detail::default_initial_vector(a.rows()), options);
}

inline EigenResult inverse_power_eigenpair(const DenseMatrix& a,
                                           const Vector& initial_vector,
                                           const EigenOptions& options) {
    EigenMetrics metrics;
    auto setup = eigen_detail::validate_setup(a, &initial_vector, options, "inverse_power", metrics);
    if (!setup.success()) {
        return setup;
    }
    if (eigen_detail::zero_vector(initial_vector)) {
        return eigen_detail::make_result(EigenStatus::invalid_input, "inverse_power", "zero_initial_vector",
                                         "invalid_input", "validation", "Initial eigenvector must be nonzero", metrics);
    }

    DenseMatrix shifted = a;
    for (std::size_t i = 0; i < shifted.rows(); ++i) {
        shifted(i, i) -= options.shift;
    }

    Vector x = initial_vector;
    const double initial_scale = eigen_detail::vector_scale(x);
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] /= initial_scale;
    }

    double previous_eigenvalue = 0.0;
    double eigenvalue = options.shift;
    for (std::size_t iteration = 1; iteration <= options.max_iterations; ++iteration) {
        auto linear = solve_dense_partial_pivot(shifted, x, options.linear_options);
        ++metrics.linear_solves;
        metrics.operation_count += linear.metrics.operation_count;
        if (!linear.success()) {
            const auto status = linear.status == SolverStatus::singular ? EigenStatus::singular_shift : EigenStatus::breakdown;
            return eigen_detail::make_result(status, "inverse_power", "linear_solve_failed",
                                             linear.status == SolverStatus::singular ? "singular_shift" : "linear_solve_failed",
                                             "linear_solve",
                                             std::string("Inverse power shifted solve failed: ") + to_string(linear.status),
                                             metrics);
        }
        Vector y = linear.solution;
        const double y_scale = eigen_detail::vector_scale(y);
        if (y_scale == 0.0 || !std::isfinite(y_scale)) {
            return eigen_detail::make_result(EigenStatus::breakdown, "inverse_power", "zero_inverse_vector",
                                             "zero_or_non_finite_vector", "normalization",
                                             "Inverse power iteration produced a zero or non-finite vector", metrics);
        }
        const std::size_t pivot = eigen_detail::max_abs_index(y);
        const double signed_scale = y[pivot];
        for (std::size_t i = 0; i < y.size(); ++i) {
            y[i] /= signed_scale;
        }
        x = std::move(y);

        try {
            eigenvalue = eigen_detail::rayleigh_quotient(a, x, metrics.operation_count);
            metrics.residual_norm = eigen_detail::eigen_residual(a, x, eigenvalue, metrics.operation_count);
        } catch (const std::exception& e) {
            return eigen_detail::make_result(EigenStatus::breakdown, "inverse_power", "rayleigh_breakdown",
                                             "non_finite_intermediate", "rayleigh", e.what(), metrics);
        }
        metrics.iterations = iteration;
        metrics.shift = options.shift;
        metrics.eigenvalue_delta = std::abs(eigenvalue - previous_eigenvalue);
        metrics.eigenvector_scale = eigen_detail::vector_scale(x);
        metrics.residual_acceptance_tolerance =
            options.comparison_tolerance(metrics.matrix_scale * metrics.eigenvector_scale +
                                         std::abs(eigenvalue) * metrics.eigenvector_scale);
        if (metrics.residual_norm <= metrics.residual_acceptance_tolerance) {
            return eigen_detail::finish_success("inverse_power", eigenvalue, x, metrics);
        }
        previous_eigenvalue = eigenvalue;
    }

    return eigen_detail::make_result(EigenStatus::not_converged, "inverse_power", "maximum_iterations",
                                     "maximum_iterations", "convergence",
                                     "Inverse power iteration reached maximum iterations", metrics);
}

} // namespace MatCal::Linalg

#endif
