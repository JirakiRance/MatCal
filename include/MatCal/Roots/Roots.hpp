#ifndef MATCAL_ROOTS_ROOTS_HPP
#define MATCAL_ROOTS_ROOTS_HPP

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace MatCal::Roots {

enum class RootStatus {
    success,
    invalid_input,
    bracket_error,
    derivative_zero,
    denominator_zero,
    non_finite,
    not_converged
};

enum class RootReason {
    none,
    invalid_callable,
    invalid_options,
    invalid_interval,
    invalid_initial_value,
    no_sign_change,
    zero_derivative,
    zero_denominator,
    non_finite_input,
    non_finite_function_value,
    non_finite_iterate,
    maximum_iterations
};

enum class RootPhase {
    setup,
    bracket,
    iteration,
    derivative,
    denominator,
    convergence
};

struct RootOptions {
    double absolute_tolerance = 1.0e-10;
    double relative_tolerance = 1.0e-10;
    double derivative_tolerance = 1.0e-14;
    int max_iterations = 1000;
    int downhill_max_shrinks = 20;
    double divergence_growth_factor = 1.5;
    int divergence_limit = 5;

    bool valid() const noexcept {
        return std::isfinite(absolute_tolerance) &&
               std::isfinite(relative_tolerance) &&
               std::isfinite(derivative_tolerance) &&
               absolute_tolerance >= 0.0 &&
               relative_tolerance >= 0.0 &&
               derivative_tolerance >= 0.0 &&
               max_iterations >= 0 &&
               downhill_max_shrinks >= 0 &&
               std::isfinite(divergence_growth_factor) &&
               divergence_growth_factor >= 1.0 &&
               divergence_limit >= 0 &&
               (absolute_tolerance > 0.0 || relative_tolerance > 0.0);
    }
};

struct RootDiagnostic {
    RootStatus status = RootStatus::success;
    RootReason reason = RootReason::none;
    RootPhase phase = RootPhase::setup;
    int iteration = 0;
    double value = 0.0;
    std::string message;
};

struct RootMetrics {
    int iterations = 0;
    double function_value = 0.0;
    double residual = 0.0;
    double final_step = 0.0;
    double absolute_tolerance_used = 0.0;
    double relative_tolerance_used = 0.0;
    double tolerance_used = 0.0;
};

struct RootResult {
    double value = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;
    RootDiagnostic diagnostic;
    RootMetrics metrics;
    std::vector<double> series;

    bool success() const noexcept {
        return converged && diagnostic.status == RootStatus::success;
    }
};

using Function = std::function<double(double)>;

namespace detail {

inline bool finite(double value) noexcept {
    return std::isfinite(value);
}

inline double tolerance_for(double value, const RootOptions& options) noexcept {
    return std::max(options.absolute_tolerance, options.relative_tolerance * std::max(std::abs(value), 1.0));
}

inline RootResult make_failure(RootStatus status,
                               RootReason reason,
                               RootPhase phase,
                               const RootOptions& options,
                               const std::string& message,
                               int iteration = 0,
                               double value = std::numeric_limits<double>::quiet_NaN(),
                               double function_value = std::numeric_limits<double>::quiet_NaN(),
                               double step = std::numeric_limits<double>::quiet_NaN(),
                               std::vector<double> series = {}) {
    RootResult result;
    result.value = value;
    result.converged = false;
    result.diagnostic = {status, reason, phase, iteration, value, message};
    result.metrics.iterations = iteration;
    result.metrics.function_value = function_value;
    result.metrics.residual = std::abs(function_value);
    result.metrics.final_step = step;
    result.metrics.absolute_tolerance_used = options.absolute_tolerance;
    result.metrics.relative_tolerance_used = options.relative_tolerance;
    result.metrics.tolerance_used = finite(value) ? tolerance_for(value, options) : options.absolute_tolerance;
    result.series = std::move(series);
    return result;
}

inline RootResult make_success(double value,
                               double function_value,
                               double step,
                               int iteration,
                               const RootOptions& options,
                               std::vector<double> series) {
    RootResult result;
    result.value = value;
    result.converged = true;
    result.diagnostic = {RootStatus::success, RootReason::none, RootPhase::convergence, iteration, value, "converged"};
    result.metrics.iterations = iteration;
    result.metrics.function_value = function_value;
    result.metrics.residual = std::abs(function_value);
    result.metrics.final_step = std::abs(step);
    result.metrics.absolute_tolerance_used = options.absolute_tolerance;
    result.metrics.relative_tolerance_used = options.relative_tolerance;
    result.metrics.tolerance_used = tolerance_for(value, options);
    result.series = std::move(series);
    return result;
}

inline RootResult validate_common(const Function& f, const RootOptions& options) {
    if (!f) {
        return make_failure(RootStatus::invalid_input, RootReason::invalid_callable, RootPhase::setup,
                            options, "root function is empty");
    }
    if (!options.valid()) {
        return make_failure(RootStatus::invalid_input, RootReason::invalid_options, RootPhase::setup,
                            options, "root options are invalid");
    }
    return {};
}

inline bool common_failed(const RootResult& result) noexcept {
    return result.diagnostic.status != RootStatus::success || result.diagnostic.reason != RootReason::none || !result.series.empty();
}

inline bool converged(double value, double function_value, double step, const RootOptions& options) noexcept {
    const double tolerance = tolerance_for(value, options);
    return std::abs(function_value) <= tolerance || std::abs(step) <= tolerance;
}

inline double checked_eval(const Function& f, double x, RootResult& failure, const RootOptions& options,
                           RootPhase phase, int iteration, std::vector<double> series = {}) {
    if (!finite(x)) {
        failure = make_failure(RootStatus::non_finite, RootReason::non_finite_input, phase, options,
                               "root evaluation point is not finite", iteration, x, std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::quiet_NaN(), std::move(series));
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double fx = f(x);
    if (!finite(fx)) {
        failure = make_failure(RootStatus::non_finite, RootReason::non_finite_function_value, phase, options,
                               "root function returned a non-finite value", iteration, x, fx,
                               std::numeric_limits<double>::quiet_NaN(), std::move(series));
    }
    return fx;
}

} // namespace detail

inline RootResult solve_bisection(const Function& f, double a, double b, RootOptions options = {}) {
    RootResult common = detail::validate_common(f, options);
    if (detail::common_failed(common)) {
        return common;
    }
    if (!detail::finite(a) || !detail::finite(b)) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::non_finite_input, RootPhase::setup,
                                    options, "bisection interval endpoints must be finite");
    }
    if (!(a < b)) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::invalid_interval, RootPhase::setup,
                                    options, "bisection interval must satisfy a < b");
    }

    RootResult failure;
    double fa = detail::checked_eval(f, a, failure, options, RootPhase::bracket, 0);
    if (detail::common_failed(failure)) return failure;
    double fb = detail::checked_eval(f, b, failure, options, RootPhase::bracket, 0);
    if (detail::common_failed(failure)) return failure;
    std::vector<double> series{a, b};

    if (fa == 0.0) return detail::make_success(a, fa, 0.0, 0, options, series);
    if (fb == 0.0) return detail::make_success(b, fb, 0.0, 0, options, series);
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        return detail::make_failure(RootStatus::bracket_error, RootReason::no_sign_change, RootPhase::bracket,
                                    options, "bisection interval endpoints do not bracket a sign change",
                                    0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                                    b - a, series);
    }

    double midpoint = a;
    double fmid = fa;
    double step = b - a;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        midpoint = a + (b - a) / 2.0;
        fmid = detail::checked_eval(f, midpoint, failure, options, RootPhase::iteration, iteration, series);
        if (detail::common_failed(failure)) return failure;
        series.push_back(midpoint);
        step = (b - a) / 2.0;
        if (fmid == 0.0 || std::abs(step) <= detail::tolerance_for(midpoint, options)) {
            return detail::make_success(midpoint, fmid, step, iteration, options, series);
        }
        if ((fa > 0.0 && fmid < 0.0) || (fa < 0.0 && fmid > 0.0)) {
            b = midpoint;
            fb = fmid;
        } else {
            a = midpoint;
            fa = fmid;
        }
        (void)fb;
    }
    return detail::make_failure(RootStatus::not_converged, RootReason::maximum_iterations, RootPhase::iteration,
                                options, "bisection reached maximum iterations", options.max_iterations,
                                midpoint, fmid, step, series);
}

inline RootResult solve_picard(const Function& phi, double x0, RootOptions options = {}) {
    RootResult common = detail::validate_common(phi, options);
    if (detail::common_failed(common)) return common;
    if (!detail::finite(x0)) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::invalid_initial_value, RootPhase::setup,
                                    options, "Picard initial value must be finite");
    }

    std::vector<double> series{x0};
    double current = x0;
    double last_step = std::numeric_limits<double>::infinity();
    int growth_count = 0;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        RootResult failure;
        const double next = detail::checked_eval(phi, current, failure, options, RootPhase::iteration, iteration, series);
        if (detail::common_failed(failure)) return failure;
        const double step = next - current;
        series.push_back(next);
        if (!detail::finite(step)) {
            return detail::make_failure(RootStatus::non_finite, RootReason::non_finite_iterate, RootPhase::iteration,
                                        options, "Picard step became non-finite", iteration, next, step, step, series);
        }
        if (std::abs(step) <= detail::tolerance_for(next, options)) {
            return detail::make_success(next, step, step, iteration, options, series);
        }
        if (std::abs(step) > last_step * options.divergence_growth_factor) {
            ++growth_count;
            if (growth_count >= options.divergence_limit && options.divergence_limit > 0) {
                return detail::make_failure(RootStatus::not_converged, RootReason::maximum_iterations, RootPhase::iteration,
                                            options, "Picard iteration appears divergent", iteration, next, step, step, series);
            }
        } else {
            growth_count = 0;
        }
        last_step = std::abs(step);
        current = next;
    }
    const double residual = series.size() >= 2 ? series.back() - series[series.size() - 2] : 0.0;
    return detail::make_failure(RootStatus::not_converged, RootReason::maximum_iterations, RootPhase::iteration,
                                options, "Picard reached maximum iterations", options.max_iterations,
                                current, residual, residual, series);
}

inline RootResult solve_picard_aitken(const Function& phi, double x0, RootOptions options = {}) {
    if (!phi) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::invalid_callable, RootPhase::setup,
                                    options, "Picard-Aitken function is empty");
    }
    Function aitken = [phi, options](double x) {
        const double p1 = phi(x);
        const double p2 = phi(p1);
        const double denominator = p2 - 2.0 * p1 + x;
        if (!std::isfinite(p1) || !std::isfinite(p2) || !std::isfinite(denominator)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (std::abs(denominator) <= options.derivative_tolerance) {
            return p1;
        }
        return (x * p2 - p1 * p1) / denominator;
    };
    return solve_picard(aitken, x0, options);
}

inline RootResult solve_newton(const Function& f, const Function& derivative, double x0, RootOptions options = {}) {
    RootResult common = detail::validate_common(f, options);
    if (detail::common_failed(common)) return common;
    if (!derivative) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::invalid_callable, RootPhase::setup,
                                    options, "Newton derivative function is empty");
    }
    if (!detail::finite(x0)) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::invalid_initial_value, RootPhase::setup,
                                    options, "Newton initial value must be finite");
    }

    std::vector<double> series{x0};
    double current = x0;
    double f_current = 0.0;
    double last_step = std::numeric_limits<double>::infinity();
    int growth_count = 0;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        RootResult failure;
        f_current = detail::checked_eval(f, current, failure, options, RootPhase::iteration, iteration, series);
        if (detail::common_failed(failure)) return failure;
        if (std::abs(f_current) <= detail::tolerance_for(current, options)) {
            return detail::make_success(current, f_current, 0.0, iteration - 1, options, series);
        }
        const double df_current = detail::checked_eval(derivative, current, failure, options, RootPhase::derivative, iteration, series);
        if (detail::common_failed(failure)) return failure;
        if (std::abs(df_current) <= options.derivative_tolerance) {
            return detail::make_failure(RootStatus::derivative_zero, RootReason::zero_derivative, RootPhase::derivative,
                                        options, "Newton derivative is near zero", iteration, current, f_current, 0.0, series);
        }
        const double step = -f_current / df_current;
        const double next = current + step;
        if (!detail::finite(next) || !detail::finite(step)) {
            return detail::make_failure(RootStatus::non_finite, RootReason::non_finite_iterate, RootPhase::iteration,
                                        options, "Newton iterate became non-finite", iteration, next, f_current, step, series);
        }
        RootResult next_failure;
        const double f_next = detail::checked_eval(f, next, next_failure, options, RootPhase::iteration, iteration, series);
        if (detail::common_failed(next_failure)) return next_failure;
        series.push_back(next);
        if (detail::converged(next, f_next, step, options)) {
            return detail::make_success(next, f_next, step, iteration, options, series);
        }
        if (std::abs(step) > last_step) {
            ++growth_count;
            if (growth_count >= options.divergence_limit && options.divergence_limit > 0) {
                return detail::make_failure(RootStatus::not_converged, RootReason::maximum_iterations, RootPhase::iteration,
                                            options, "Newton iteration appears divergent", iteration, next, f_next, step, series);
            }
        } else {
            growth_count = 0;
        }
        last_step = std::abs(step);
        current = next;
        f_current = f_next;
    }
    return detail::make_failure(RootStatus::not_converged, RootReason::maximum_iterations, RootPhase::iteration,
                                options, "Newton reached maximum iterations", options.max_iterations,
                                current, f_current, 0.0, series);
}

inline RootResult solve_downhill_newton(const Function& f, const Function& derivative, double x0, RootOptions options = {}) {
    RootResult common = detail::validate_common(f, options);
    if (detail::common_failed(common)) return common;
    if (!derivative) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::invalid_callable, RootPhase::setup,
                                    options, "Downhill Newton derivative function is empty");
    }
    if (!detail::finite(x0)) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::invalid_initial_value, RootPhase::setup,
                                    options, "Downhill Newton initial value must be finite");
    }

    std::vector<double> series{x0};
    double current = x0;
    RootResult failure;
    double f_current = detail::checked_eval(f, current, failure, options, RootPhase::iteration, 0, series);
    if (detail::common_failed(failure)) return failure;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        if (std::abs(f_current) <= detail::tolerance_for(current, options)) {
            return detail::make_success(current, f_current, 0.0, iteration - 1, options, series);
        }
        const double df_current = detail::checked_eval(derivative, current, failure, options, RootPhase::derivative, iteration, series);
        if (detail::common_failed(failure)) return failure;
        if (std::abs(df_current) <= options.derivative_tolerance) {
            return detail::make_failure(RootStatus::derivative_zero, RootReason::zero_derivative, RootPhase::derivative,
                                        options, "Downhill Newton derivative is near zero", iteration, current, f_current, 0.0, series);
        }
        const double newton_step = -f_current / df_current;
        double lambda = 1.0;
        double next = current + newton_step;
        double f_next = detail::checked_eval(f, next, failure, options, RootPhase::iteration, iteration, series);
        if (detail::common_failed(failure)) return failure;
        int shrink = 0;
        while (std::abs(f_next) > std::abs(f_current) && shrink < options.downhill_max_shrinks) {
            lambda *= 0.5;
            next = current + lambda * newton_step;
            f_next = detail::checked_eval(f, next, failure, options, RootPhase::iteration, iteration, series);
            if (detail::common_failed(failure)) return failure;
            ++shrink;
        }
        const double step = next - current;
        if (!detail::finite(next) || !detail::finite(step)) {
            return detail::make_failure(RootStatus::non_finite, RootReason::non_finite_iterate, RootPhase::iteration,
                                        options, "Downhill Newton iterate became non-finite", iteration, next, f_next, step, series);
        }
        series.push_back(next);
        if (detail::converged(next, f_next, step, options)) {
            return detail::make_success(next, f_next, step, iteration, options, series);
        }
        current = next;
        f_current = f_next;
    }
    return detail::make_failure(RootStatus::not_converged, RootReason::maximum_iterations, RootPhase::iteration,
                                options, "Downhill Newton reached maximum iterations", options.max_iterations,
                                current, f_current, 0.0, series);
}

inline RootResult solve_secant_two_point(const Function& f, double x0, double x1, RootOptions options = {}) {
    RootResult common = detail::validate_common(f, options);
    if (detail::common_failed(common)) return common;
    if (!detail::finite(x0) || !detail::finite(x1)) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::invalid_initial_value, RootPhase::setup,
                                    options, "secant initial values must be finite");
    }
    std::vector<double> series{x0, x1};
    RootResult failure;
    double f0 = detail::checked_eval(f, x0, failure, options, RootPhase::iteration, 0, series);
    if (detail::common_failed(failure)) return failure;
    double f1 = detail::checked_eval(f, x1, failure, options, RootPhase::iteration, 0, series);
    if (detail::common_failed(failure)) return failure;
    if (std::abs(f0) <= detail::tolerance_for(x0, options)) return detail::make_success(x0, f0, 0.0, 0, options, series);
    if (std::abs(f1) <= detail::tolerance_for(x1, options)) return detail::make_success(x1, f1, 0.0, 0, options, series);

    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        const double denominator = f1 - f0;
        if (std::abs(denominator) <= options.derivative_tolerance) {
            return detail::make_failure(RootStatus::denominator_zero, RootReason::zero_denominator, RootPhase::denominator,
                                        options, "secant denominator is near zero", iteration, x1, f1, x1 - x0, series);
        }
        const double next = x1 - f1 * (x1 - x0) / denominator;
        if (!detail::finite(next)) {
            return detail::make_failure(RootStatus::non_finite, RootReason::non_finite_iterate, RootPhase::iteration,
                                        options, "secant iterate became non-finite", iteration, next, f1, x1 - x0, series);
        }
        const double f_next = detail::checked_eval(f, next, failure, options, RootPhase::iteration, iteration, series);
        if (detail::common_failed(failure)) return failure;
        const double step = next - x1;
        series.push_back(next);
        if (detail::converged(next, f_next, step, options)) {
            return detail::make_success(next, f_next, step, iteration, options, series);
        }
        x0 = x1;
        f0 = f1;
        x1 = next;
        f1 = f_next;
    }
    return detail::make_failure(RootStatus::not_converged, RootReason::maximum_iterations, RootPhase::iteration,
                                options, "secant reached maximum iterations", options.max_iterations,
                                x1, f1, 0.0, series);
}

inline RootResult solve_secant_one_point(const Function& f, double x0, double h = 1.0e-4, RootOptions options = {}) {
    if (!detail::finite(h) || h == 0.0) {
        return detail::make_failure(RootStatus::invalid_input, RootReason::invalid_initial_value, RootPhase::setup,
                                    options, "one-point secant step h must be finite and non-zero");
    }
    return solve_secant_two_point(f, x0 - h, x0, options);
}

} // namespace MatCal::Roots

#endif
