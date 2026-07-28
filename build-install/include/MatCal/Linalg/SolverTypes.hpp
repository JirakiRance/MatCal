#ifndef MATCAL_LINALG_SOLVER_TYPES_HPP
#define MATCAL_LINALG_SOLVER_TYPES_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "MatCal/Linalg/Vector.hpp"

namespace MatCal::Linalg {

enum class SolverStatus {
    success,
    invalid_input,
    dimension_mismatch,
    non_finite_input,
    singular,
    not_positive_definite,
    breakdown,
    not_converged
};

inline const char* to_string(SolverStatus status) noexcept {
    switch (status) {
    case SolverStatus::success: return "success";
    case SolverStatus::invalid_input: return "invalid_input";
    case SolverStatus::dimension_mismatch: return "dimension_mismatch";
    case SolverStatus::non_finite_input: return "non_finite_input";
    case SolverStatus::singular: return "singular";
    case SolverStatus::not_positive_definite: return "not_positive_definite";
    case SolverStatus::breakdown: return "breakdown";
    case SolverStatus::not_converged: return "not_converged";
    }
    return "unknown";
}

struct SolverOptions {
    double absolute_tolerance = 0.0;
    double relative_tolerance = 1e-12;
    double pivot_factor = 1.0;
    std::size_t max_iterations = 1000;

    bool valid() const noexcept {
        return std::isfinite(absolute_tolerance) &&
               std::isfinite(relative_tolerance) &&
               std::isfinite(pivot_factor) &&
               absolute_tolerance >= 0.0 &&
               relative_tolerance >= 0.0 &&
               pivot_factor >= 0.0 &&
               max_iterations > 0;
    }

    double comparison_tolerance(double scale) const noexcept {
        if (!std::isfinite(scale) || scale < 0.0) {
            return std::numeric_limits<double>::max();
        }
        double relative_part = saturated_product(relative_tolerance, scale);
        return std::max(absolute_tolerance, relative_part);
    }

    double pivot_tolerance(double scale) const noexcept {
        return saturated_product(pivot_factor, comparison_tolerance(scale));
    }

private:
    static double saturated_product(double left, double right) noexcept {
        if (left == 0.0 || right == 0.0) {
            return 0.0;
        }
        double max = std::numeric_limits<double>::max();
        if (left > max / right) {
            return max;
        }
        return left * right;
    }

    static double saturated_sum(double left, double right) noexcept {
        double max = std::numeric_limits<double>::max();
        if (left > max - right) {
            return max;
        }
        return left + right;
    }
};

struct SolverDiagnostic {
    SolverStatus status = SolverStatus::success;
    std::string code;
    std::string reason;
    std::string phase;
    std::size_t row = invalid_index();
    std::size_t column = invalid_index();
    double value = 0.0;
    double scale = 0.0;
    double tolerance = 0.0;
    std::string message;

    static constexpr std::size_t invalid_index() noexcept {
        return std::numeric_limits<std::size_t>::max();
    }
};

struct SolverMetrics {
    std::size_t iterations = 0;
    std::size_t operation_count = 0;
    std::size_t factorization_operation_count = 0;
    std::size_t solve_operation_count = 0;
    double residual_norm = 0.0;
    double absolute_residual_norm = 0.0;
    double relative_residual_norm = 0.0;
    double residual_acceptance_tolerance = 0.0;
    double matrix_scale = 0.0;
    double rhs_scale = 0.0;
    double solution_scale = 0.0;
    double pivot_tolerance_used = 0.0;
    double minimum_abs_pivot = std::numeric_limits<double>::infinity();
};

struct SolverResult {
    SolverStatus status = SolverStatus::invalid_input;
    Vector solution;
    SolverMetrics metrics;
    std::vector<SolverDiagnostic> diagnostics;
    std::string method;
    std::string implementation;

    bool success() const noexcept {
        return status == SolverStatus::success;
    }
};

} // namespace MatCal::Linalg

#endif
