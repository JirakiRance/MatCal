#ifndef MATCAL_LINALG_SOLVER_TYPES_HPP
#define MATCAL_LINALG_SOLVER_TYPES_HPP

#include <cmath>
#include <cstddef>
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
    double absolute_tolerance = 1e-12;
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
        return absolute_tolerance + relative_tolerance * scale;
    }

    double pivot_tolerance(double scale) const noexcept {
        return pivot_factor * comparison_tolerance(scale);
    }
};

struct SolverDiagnostic {
    SolverStatus status = SolverStatus::success;
    std::string message;
};

struct SolverMetrics {
    std::size_t iterations = 0;
    double residual_norm = 0.0;
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
