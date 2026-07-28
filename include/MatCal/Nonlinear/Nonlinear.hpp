#ifndef MATCAL_NONLINEAR_NONLINEAR_HPP
#define MATCAL_NONLINEAR_NONLINEAR_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/SolverTypes.hpp"
#include "MatCal/Linalg/Vector.hpp"

namespace MatCal::Nonlinear {

enum class NonlinearStatus {
    success,
    invalid_input,
    dimension_mismatch,
    non_finite_input,
    singular_jacobian,
    breakdown,
    not_converged
};

enum class NonlinearReason {
    none,
    invalid_callable,
    invalid_options,
    dimension_mismatch,
    non_finite_initial_guess,
    non_finite_residual,
    non_finite_jacobian,
    non_finite_step,
    singular_jacobian,
    linear_solve_failed,
    maximum_iterations
};

enum class NonlinearPhase {
    setup,
    residual,
    jacobian,
    linear_solve,
    update,
    convergence
};

struct NonlinearOptions {
    double absolute_tolerance = 1.0e-10;
    double relative_tolerance = 1.0e-10;
    double finite_difference_step = 1.0e-6;
    std::size_t max_iterations = 50;
    MatCal::Linalg::SolverOptions linear_options{};

    bool valid() const noexcept {
        return std::isfinite(absolute_tolerance) &&
               std::isfinite(relative_tolerance) &&
               std::isfinite(finite_difference_step) &&
               absolute_tolerance >= 0.0 &&
               relative_tolerance >= 0.0 &&
               finite_difference_step > 0.0 &&
               max_iterations > 0 &&
               linear_options.valid();
    }
};

struct NonlinearDiagnostic {
    NonlinearStatus status = NonlinearStatus::success;
    NonlinearReason reason = NonlinearReason::none;
    NonlinearPhase phase = NonlinearPhase::setup;
    std::size_t row = invalid_index();
    std::size_t column = invalid_index();
    double value = 0.0;
    std::string message;

    static constexpr std::size_t invalid_index() noexcept {
        return std::numeric_limits<std::size_t>::max();
    }
};

struct NonlinearMetrics {
    std::size_t iterations = 0;
    std::size_t function_evaluations = 0;
    std::size_t jacobian_evaluations = 0;
    std::size_t linear_solves = 0;
    double residual_norm = std::numeric_limits<double>::infinity();
    double step_norm = std::numeric_limits<double>::infinity();
    double absolute_tolerance = 0.0;
    double relative_tolerance = 0.0;
};

struct NonlinearResult {
    NonlinearStatus status = NonlinearStatus::invalid_input;
    std::vector<double> solution;
    bool converged = false;
    NonlinearDiagnostic diagnostic;
    NonlinearMetrics metrics;

    bool success() const noexcept {
        return status == NonlinearStatus::success && converged;
    }
};

using ResidualFunction = std::function<std::vector<double>(const std::vector<double>&)>;
using JacobianFunction = std::function<MatCal::Linalg::DenseMatrix(const std::vector<double>&)>;

namespace detail {

inline bool finite_vector(const std::vector<double>& values) noexcept {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

inline double norm_inf(const std::vector<double>& values) {
    double result = 0.0;
    for (double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

inline NonlinearDiagnostic diagnostic(NonlinearStatus status,
                                      NonlinearReason reason,
                                      NonlinearPhase phase,
                                      std::string message,
                                      std::size_t row = NonlinearDiagnostic::invalid_index(),
                                      std::size_t column = NonlinearDiagnostic::invalid_index(),
                                      double value = 0.0) {
    return {status, reason, phase, row, column, value, std::move(message)};
}

inline NonlinearResult failure(NonlinearStatus status,
                               NonlinearReason reason,
                               NonlinearPhase phase,
                               std::string message,
                               NonlinearMetrics metrics,
                               std::size_t row = NonlinearDiagnostic::invalid_index(),
                               std::size_t column = NonlinearDiagnostic::invalid_index(),
                               double value = 0.0) {
    NonlinearResult result;
    result.status = status;
    result.converged = false;
    result.metrics = metrics;
    result.diagnostic = diagnostic(status, reason, phase, std::move(message), row, column, value);
    return result;
}

inline bool residual_accepts(double residual_norm, double initial_residual_norm, const NonlinearOptions& options) {
    const double scale = std::max(initial_residual_norm, 1.0);
    return residual_norm <= std::max(options.absolute_tolerance, options.relative_tolerance * scale);
}

inline bool step_accepts(double step_norm, double x_norm, const NonlinearOptions& options) {
    const double scale = std::max(x_norm, 1.0);
    return step_norm <= std::max(options.absolute_tolerance, options.relative_tolerance * scale);
}

inline MatCal::Linalg::Vector to_vector(const std::vector<double>& values) {
    MatCal::Linalg::Vector result(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        result[i] = values[i];
    }
    return result;
}

inline std::vector<double> from_vector(const MatCal::Linalg::Vector& values) {
    std::vector<double> result(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        result[i] = values[i];
    }
    return result;
}

} // namespace detail

inline NonlinearResult solve_newton_system(const ResidualFunction& residual,
                                           const std::vector<double>& initial_guess,
                                           const JacobianFunction& jacobian,
                                           NonlinearOptions options = {}) {
    NonlinearMetrics metrics;
    metrics.absolute_tolerance = options.absolute_tolerance;
    metrics.relative_tolerance = options.relative_tolerance;

    if (!residual) {
        return detail::failure(NonlinearStatus::invalid_input, NonlinearReason::invalid_callable,
                               NonlinearPhase::setup, "nonlinear residual callable is empty", metrics);
    }
    if (!jacobian) {
        return detail::failure(NonlinearStatus::invalid_input, NonlinearReason::invalid_callable,
                               NonlinearPhase::setup, "nonlinear jacobian callable is empty", metrics);
    }
    if (!options.valid()) {
        return detail::failure(NonlinearStatus::invalid_input, NonlinearReason::invalid_options,
                               NonlinearPhase::setup, "nonlinear options are invalid", metrics);
    }
    if (initial_guess.empty()) {
        return detail::failure(NonlinearStatus::invalid_input, NonlinearReason::dimension_mismatch,
                               NonlinearPhase::setup, "nonlinear initial guess must not be empty", metrics);
    }
    if (!detail::finite_vector(initial_guess)) {
        return detail::failure(NonlinearStatus::non_finite_input, NonlinearReason::non_finite_initial_guess,
                               NonlinearPhase::setup, "nonlinear initial guess must be finite", metrics);
    }

    std::vector<double> x = initial_guess;
    const std::size_t n = x.size();

    std::vector<double> f;
    try {
        f = residual(x);
        ++metrics.function_evaluations;
    } catch (const std::exception& e) {
        return detail::failure(NonlinearStatus::breakdown, NonlinearReason::non_finite_residual,
                               NonlinearPhase::residual, e.what(), metrics);
    }
    if (f.size() != n) {
        return detail::failure(NonlinearStatus::dimension_mismatch, NonlinearReason::dimension_mismatch,
                               NonlinearPhase::residual, "nonlinear residual dimension mismatch", metrics);
    }
    if (!detail::finite_vector(f)) {
        return detail::failure(NonlinearStatus::breakdown, NonlinearReason::non_finite_residual,
                               NonlinearPhase::residual, "nonlinear residual produced a non-finite value", metrics);
    }
    const double initial_residual_norm = detail::norm_inf(f);
    metrics.residual_norm = initial_residual_norm;
    if (detail::residual_accepts(metrics.residual_norm, initial_residual_norm, options)) {
        return {NonlinearStatus::success, x, true,
                detail::diagnostic(NonlinearStatus::success, NonlinearReason::none,
                                   NonlinearPhase::convergence, "nonlinear solve converged at initial guess"),
                metrics};
    }

    for (std::size_t iteration = 1; iteration <= options.max_iterations; ++iteration) {
        metrics.iterations = iteration;
        MatCal::Linalg::DenseMatrix j;
        try {
            j = jacobian(x);
            ++metrics.jacobian_evaluations;
        } catch (const std::exception& e) {
            return detail::failure(NonlinearStatus::breakdown, NonlinearReason::non_finite_jacobian,
                                   NonlinearPhase::jacobian, e.what(), metrics);
        }
        if (j.rows() != n || j.cols() != n) {
            return detail::failure(NonlinearStatus::dimension_mismatch, NonlinearReason::dimension_mismatch,
                                   NonlinearPhase::jacobian, "nonlinear jacobian dimension mismatch", metrics);
        }
        if (!j.all_finite()) {
            return detail::failure(NonlinearStatus::breakdown, NonlinearReason::non_finite_jacobian,
                                   NonlinearPhase::jacobian, "nonlinear jacobian produced a non-finite value", metrics);
        }

        std::vector<double> rhs_values(n);
        for (std::size_t i = 0; i < n; ++i) {
            rhs_values[i] = -f[i];
        }

        auto linear_result = MatCal::Linalg::solve_dense_partial_pivot(j, detail::to_vector(rhs_values), options.linear_options);
        ++metrics.linear_solves;
        if (!linear_result.success()) {
            const auto status = linear_result.status == MatCal::Linalg::SolverStatus::singular
                                    ? NonlinearStatus::singular_jacobian
                                    : NonlinearStatus::breakdown;
            const auto reason = linear_result.status == MatCal::Linalg::SolverStatus::singular
                                    ? NonlinearReason::singular_jacobian
                                    : NonlinearReason::linear_solve_failed;
            return detail::failure(status, reason, NonlinearPhase::linear_solve,
                                   std::string("nonlinear Newton linear solve failed: ") + MatCal::Linalg::to_string(linear_result.status),
                                   metrics);
        }

        std::vector<double> step = detail::from_vector(linear_result.solution);
        if (!detail::finite_vector(step)) {
            return detail::failure(NonlinearStatus::breakdown, NonlinearReason::non_finite_step,
                                   NonlinearPhase::linear_solve, "nonlinear Newton step is non-finite", metrics);
        }
        metrics.step_norm = detail::norm_inf(step);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += step[i];
            if (!std::isfinite(x[i])) {
                return detail::failure(NonlinearStatus::breakdown, NonlinearReason::non_finite_step,
                                       NonlinearPhase::update, "nonlinear Newton update produced a non-finite value",
                                       metrics, i, NonlinearDiagnostic::invalid_index(), x[i]);
            }
        }

        try {
            f = residual(x);
            ++metrics.function_evaluations;
        } catch (const std::exception& e) {
            return detail::failure(NonlinearStatus::breakdown, NonlinearReason::non_finite_residual,
                                   NonlinearPhase::residual, e.what(), metrics);
        }
        if (f.size() != n) {
            return detail::failure(NonlinearStatus::dimension_mismatch, NonlinearReason::dimension_mismatch,
                                   NonlinearPhase::residual, "nonlinear residual dimension mismatch after update", metrics);
        }
        if (!detail::finite_vector(f)) {
            return detail::failure(NonlinearStatus::breakdown, NonlinearReason::non_finite_residual,
                                   NonlinearPhase::residual, "nonlinear residual produced a non-finite value after update",
                                   metrics);
        }
        metrics.residual_norm = detail::norm_inf(f);

        if (detail::residual_accepts(metrics.residual_norm, initial_residual_norm, options) ||
            (detail::step_accepts(metrics.step_norm, detail::norm_inf(x), options) &&
             detail::residual_accepts(metrics.residual_norm, initial_residual_norm, options))) {
            return {NonlinearStatus::success, x, true,
                    detail::diagnostic(NonlinearStatus::success, NonlinearReason::none,
                                       NonlinearPhase::convergence, "nonlinear solve converged"),
                    metrics};
        }
    }

    return detail::failure(NonlinearStatus::not_converged, NonlinearReason::maximum_iterations,
                           NonlinearPhase::convergence, "nonlinear solve reached maximum iterations", metrics);
}

inline MatCal::Linalg::DenseMatrix finite_difference_jacobian(const ResidualFunction& residual,
                                                             const std::vector<double>& x,
                                                             const std::vector<double>& f0,
                                                             const NonlinearOptions& options,
                                                             NonlinearMetrics& metrics) {
    const std::size_t n = x.size();
    MatCal::Linalg::DenseMatrix j(n, n);
    for (std::size_t column = 0; column < n; ++column) {
        std::vector<double> shifted = x;
        shifted[column] += options.finite_difference_step;
        if (!std::isfinite(shifted[column])) {
            throw std::overflow_error("finite-difference jacobian point is non-finite");
        }
        std::vector<double> f1 = residual(shifted);
        ++metrics.function_evaluations;
        if (f1.size() != n) {
            throw std::invalid_argument("finite-difference residual dimension mismatch");
        }
        if (!detail::finite_vector(f1)) {
            throw std::runtime_error("finite-difference residual produced non-finite value");
        }
        for (std::size_t row = 0; row < n; ++row) {
            j(row, column) = (f1[row] - f0[row]) / options.finite_difference_step;
        }
    }
    return j;
}

inline NonlinearResult solve_newton_system_finite_difference(const ResidualFunction& residual,
                                                            const std::vector<double>& initial_guess,
                                                            NonlinearOptions options = {}) {
    NonlinearMetrics shared_metrics;
    auto jacobian = [&](const std::vector<double>& x) {
        std::vector<double> f0 = residual(x);
        ++shared_metrics.function_evaluations;
        return finite_difference_jacobian(residual, x, f0, options, shared_metrics);
    };

    auto result = solve_newton_system(residual, initial_guess, jacobian, options);
    result.metrics.function_evaluations += shared_metrics.function_evaluations;
    return result;
}

} // namespace MatCal::Nonlinear

#endif
