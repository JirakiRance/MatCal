#ifndef MATCAL_CALCULUS_CALCULUS_HPP
#define MATCAL_CALCULUS_CALCULUS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace MatCal::Calculus {

enum class CalculusStatus {
    success,
    invalid_input,
    non_finite,
    not_converged
};

enum class CalculusReason {
    none,
    invalid_callable,
    invalid_step,
    invalid_interval,
    invalid_order,
    invalid_segments,
    invalid_options,
    non_finite_input,
    non_finite_function_value,
    non_finite_result,
    maximum_iterations
};

enum class CalculusPhase {
    setup,
    evaluation,
    quadrature,
    refinement,
    convergence
};

struct CalculusDiagnostic {
    CalculusStatus status = CalculusStatus::success;
    CalculusReason reason = CalculusReason::none;
    CalculusPhase phase = CalculusPhase::setup;
    int iteration = 0;
    double value = 0.0;
    std::string message;
};

struct DerivativeResult {
    double value = std::numeric_limits<double>::quiet_NaN();
    double step = 0.0;
    int function_evaluations = 0;
    CalculusDiagnostic diagnostic;

    bool success() const noexcept {
        return diagnostic.status == CalculusStatus::success;
    }
};

struct IntegrationOptions {
    double tolerance = 1.0e-10;
    int max_iterations = 20;
    int order = 4;
    int segments = 10;

    bool valid() const noexcept {
        return std::isfinite(tolerance) &&
               tolerance >= 0.0 &&
               max_iterations >= 0 &&
               order >= 1 &&
               order <= 7 &&
               segments > 0;
    }
};

struct IntegrationMetrics {
    int iterations = 0;
    int refinements = 0;
    int function_evaluations = 0;
    double tolerance = 0.0;
    double estimated_error = std::numeric_limits<double>::quiet_NaN();
};

struct IntegrationResult {
    double value = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;
    CalculusDiagnostic diagnostic;
    IntegrationMetrics metrics;
    std::vector<std::vector<double>> table;

    bool success() const noexcept {
        return converged && diagnostic.status == CalculusStatus::success;
    }
};

using Function = std::function<double(double)>;

namespace detail {

inline constexpr double nan() noexcept {
    return std::numeric_limits<double>::quiet_NaN();
}

inline bool finite(double value) noexcept {
    return std::isfinite(value);
}

inline CalculusDiagnostic success_diag() {
    return {CalculusStatus::success, CalculusReason::none, CalculusPhase::convergence, 0, 0.0, "success"};
}

inline CalculusDiagnostic failure_diag(CalculusStatus status,
                                       CalculusReason reason,
                                       CalculusPhase phase,
                                       const std::string& message,
                                       int iteration = 0,
                                       double value = nan()) {
    return {status, reason, phase, iteration, value, message};
}

inline double checked_eval(const Function& f,
                           double x,
                           CalculusDiagnostic& failure,
                           int& evaluations,
                           CalculusPhase phase = CalculusPhase::evaluation,
                           int iteration = 0) {
    if (!finite(x)) {
        failure = failure_diag(CalculusStatus::invalid_input, CalculusReason::non_finite_input,
                               phase, "evaluation point is not finite", iteration, x);
        return nan();
    }
    const double y = f(x);
    ++evaluations;
    if (!finite(y)) {
        failure = failure_diag(CalculusStatus::non_finite, CalculusReason::non_finite_function_value,
                               phase, "function returned a non-finite value", iteration, y);
    }
    return y;
}

inline DerivativeResult derivative_failure(CalculusStatus status,
                                           CalculusReason reason,
                                           const std::string& message,
                                           double step,
                                           int evaluations = 0) {
    DerivativeResult result;
    result.step = step;
    result.function_evaluations = evaluations;
    result.diagnostic = failure_diag(status, reason, CalculusPhase::setup, message);
    return result;
}

inline IntegrationResult integration_failure(CalculusStatus status,
                                             CalculusReason reason,
                                             CalculusPhase phase,
                                             const std::string& message,
                                             const IntegrationOptions& options,
                                             int evaluations = 0,
                                             int iteration = 0,
                                             double value = nan(),
                                             std::vector<std::vector<double>> table = {}) {
    IntegrationResult result;
    result.value = value;
    result.converged = false;
    result.diagnostic = failure_diag(status, reason, phase, message, iteration, value);
    result.metrics.iterations = iteration;
    result.metrics.refinements = iteration;
    result.metrics.function_evaluations = evaluations;
    result.metrics.tolerance = options.tolerance;
    result.table = std::move(table);
    return result;
}

inline IntegrationResult integration_success(double value,
                                             const IntegrationOptions& options,
                                             int evaluations,
                                             int iteration = 0,
                                             double estimated_error = nan(),
                                             std::vector<std::vector<double>> table = {}) {
    IntegrationResult result;
    result.value = value;
    result.converged = true;
    result.diagnostic = success_diag();
    result.diagnostic.iteration = iteration;
    result.diagnostic.value = value;
    result.metrics.iterations = iteration;
    result.metrics.refinements = iteration;
    result.metrics.function_evaluations = evaluations;
    result.metrics.tolerance = options.tolerance;
    result.metrics.estimated_error = estimated_error;
    result.table = std::move(table);
    return result;
}

inline double oriented_sign(double a, double b) noexcept {
    return a <= b ? 1.0 : -1.0;
}

inline std::pair<double, double> ordered_interval(double a, double b) noexcept {
    return a <= b ? std::make_pair(a, b) : std::make_pair(b, a);
}

} // namespace detail

inline const std::vector<std::vector<double>>& cotes_coefficients() {
    static const std::vector<std::vector<double>> coefficients = {
        {1.0/2, 1.0/2},
        {1.0/6, 4.0/6, 1.0/6},
        {1.0/8, 3.0/8, 3.0/8, 1.0/8},
        {7.0/90, 16.0/45, 2.0/15, 16.0/45, 7.0/90},
        {19.0/288, 25.0/96, 25.0/144, 25.0/144, 25.0/96, 19.0/288},
        {41.0/840, 9.0/35, 9.0/280, 34.0/105, 9.0/280, 9.0/35, 41.0/840},
        {751.0/17280, 3577.0/17280, 1323.0/17280, 2989.0/17280,
         2989.0/17280, 1323.0/17280, 3577.0/17280, 751.0/17280}
    };
    return coefficients;
}

inline DerivativeResult forward_difference(const Function& f, double x, double step = 1.0e-6) {
    if (!f) return detail::derivative_failure(CalculusStatus::invalid_input, CalculusReason::invalid_callable, "derivative function is empty", step);
    if (!detail::finite(x)) return detail::derivative_failure(CalculusStatus::invalid_input, CalculusReason::non_finite_input, "derivative point is not finite", step);
    if (!detail::finite(step) || step <= 0.0) return detail::derivative_failure(CalculusStatus::invalid_input, CalculusReason::invalid_step, "derivative step must be finite and positive", step);
    if (!detail::finite(x + step)) return detail::derivative_failure(CalculusStatus::invalid_input, CalculusReason::non_finite_input, "x + step is not finite", step);

    int evaluations = 0;
    CalculusDiagnostic failure;
    const double fxh = detail::checked_eval(f, x + step, failure, evaluations);
    if (failure.status != CalculusStatus::success) return {detail::nan(), step, evaluations, failure};
    const double fx = detail::checked_eval(f, x, failure, evaluations);
    if (failure.status != CalculusStatus::success) return {detail::nan(), step, evaluations, failure};
    const double value = (fxh - fx) / step;
    if (!detail::finite(value)) {
        return {detail::nan(), step, evaluations,
                detail::failure_diag(CalculusStatus::non_finite, CalculusReason::non_finite_result,
                                     CalculusPhase::evaluation, "forward difference produced a non-finite value")};
    }
    return {value, step, evaluations, detail::success_diag()};
}

inline DerivativeResult central_difference(const Function& f, double x, double step = 1.0e-6) {
    if (!f) return detail::derivative_failure(CalculusStatus::invalid_input, CalculusReason::invalid_callable, "derivative function is empty", step);
    if (!detail::finite(x)) return detail::derivative_failure(CalculusStatus::invalid_input, CalculusReason::non_finite_input, "derivative point is not finite", step);
    if (!detail::finite(step) || step <= 0.0) return detail::derivative_failure(CalculusStatus::invalid_input, CalculusReason::invalid_step, "derivative step must be finite and positive", step);
    const double half = 0.5 * step;
    if (!detail::finite(x + half) || !detail::finite(x - half)) {
        return detail::derivative_failure(CalculusStatus::invalid_input, CalculusReason::non_finite_input, "x +/- step/2 is not finite", step);
    }

    int evaluations = 0;
    CalculusDiagnostic failure;
    const double right = detail::checked_eval(f, x + half, failure, evaluations);
    if (failure.status != CalculusStatus::success) return {detail::nan(), step, evaluations, failure};
    const double left = detail::checked_eval(f, x - half, failure, evaluations);
    if (failure.status != CalculusStatus::success) return {detail::nan(), step, evaluations, failure};
    const double value = (right - left) / step;
    if (!detail::finite(value)) {
        return {detail::nan(), step, evaluations,
                detail::failure_diag(CalculusStatus::non_finite, CalculusReason::non_finite_result,
                                     CalculusPhase::evaluation, "central difference produced a non-finite value")};
    }
    return {value, step, evaluations, detail::success_diag()};
}

inline IntegrationResult integrate_instant(const Function& f, double a, double b, double step = 1.0e-6) {
    IntegrationOptions options;
    options.tolerance = step;
    if (!f) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::invalid_callable, CalculusPhase::setup, "integration function is empty", options);
    if (!detail::finite(a) || !detail::finite(b) || !detail::finite(step)) {
        return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::non_finite_input, CalculusPhase::setup, "integration interval and step must be finite", options);
    }
    if (step <= 0.0) {
        return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::invalid_step, CalculusPhase::setup, "instant integration step must be positive", options);
    }
    if (a == b) {
        return detail::integration_success(0.0, options, 0);
    }

    const double sign = detail::oriented_sign(a, b);
    auto [left, right] = detail::ordered_interval(a, b);
    double current = left;
    double total = 0.0;
    int evaluations = 0;
    while (current < right) {
        const double width = std::min(step, right - current);
        CalculusDiagnostic failure;
        const double value = detail::checked_eval(f, current, failure, evaluations, CalculusPhase::quadrature);
        if (failure.status != CalculusStatus::success) {
            return detail::integration_failure(failure.status, failure.reason, failure.phase, failure.message, options, evaluations);
        }
        total += value * width;
        if (!detail::finite(total)) {
            return detail::integration_failure(CalculusStatus::non_finite, CalculusReason::non_finite_result, CalculusPhase::quadrature, "instant integration produced a non-finite result", options, evaluations);
        }
        current += width;
    }
    return detail::integration_success(sign * total, options, evaluations);
}

inline IntegrationResult integrate_newton_cotes(const Function& f, double a, double b, int order = 4) {
    IntegrationOptions options;
    options.order = order;
    if (!f) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::invalid_callable, CalculusPhase::setup, "integration function is empty", options);
    if (!detail::finite(a) || !detail::finite(b)) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::non_finite_input, CalculusPhase::setup, "integration interval must be finite", options);
    if (order < 1 || order > 7) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::invalid_order, CalculusPhase::setup, "Newton-Cotes order must be between 1 and 7", options);
    if (a == b) return detail::integration_success(0.0, options, 0);

    const double sign = detail::oriented_sign(a, b);
    auto [left, right] = detail::ordered_interval(a, b);
    const auto& coeffs = cotes_coefficients()[static_cast<std::size_t>(order - 1)];
    const double h = (right - left) / static_cast<double>(order);
    double sum = 0.0;
    int evaluations = 0;
    for (int i = 0; i <= order; ++i) {
        CalculusDiagnostic failure;
        const double y = detail::checked_eval(f, left + i * h, failure, evaluations, CalculusPhase::quadrature);
        if (failure.status != CalculusStatus::success) return detail::integration_failure(failure.status, failure.reason, failure.phase, failure.message, options, evaluations);
        sum += coeffs[static_cast<std::size_t>(i)] * y;
    }
    const double result = sign * sum * (right - left);
    if (!detail::finite(result)) return detail::integration_failure(CalculusStatus::non_finite, CalculusReason::non_finite_result, CalculusPhase::quadrature, "Newton-Cotes produced a non-finite result", options, evaluations);
    return detail::integration_success(result, options, evaluations);
}

inline IntegrationResult integrate_composite_newton_cotes(const Function& f, double a, double b, int segments = 10, int order = 4) {
    IntegrationOptions options;
    options.order = order;
    options.segments = segments;
    if (!f) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::invalid_callable, CalculusPhase::setup, "integration function is empty", options);
    if (!detail::finite(a) || !detail::finite(b)) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::non_finite_input, CalculusPhase::setup, "integration interval must be finite", options);
    if (segments <= 0) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::invalid_segments, CalculusPhase::setup, "segments must be positive", options);
    if (order < 1 || order > 7) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::invalid_order, CalculusPhase::setup, "Newton-Cotes order must be between 1 and 7", options);
    if (a == b) return detail::integration_success(0.0, options, 0);

    const double sign = detail::oriented_sign(a, b);
    auto [left, right] = detail::ordered_interval(a, b);
    const auto& coeffs = cotes_coefficients()[static_cast<std::size_t>(order - 1)];
    const double segment_width = (right - left) / static_cast<double>(segments);
    const double h = segment_width / static_cast<double>(order);
    double total = 0.0;
    int evaluations = 0;
    for (int seg = 0; seg < segments; ++seg) {
        const double seg_a = left + seg * segment_width;
        double sum = 0.0;
        for (int i = 0; i <= order; ++i) {
            CalculusDiagnostic failure;
            const double y = detail::checked_eval(f, seg_a + i * h, failure, evaluations, CalculusPhase::quadrature);
            if (failure.status != CalculusStatus::success) return detail::integration_failure(failure.status, failure.reason, failure.phase, failure.message, options, evaluations);
            sum += coeffs[static_cast<std::size_t>(i)] * y;
        }
        total += sum;
    }
    const double result = sign * total * segment_width;
    if (!detail::finite(result)) return detail::integration_failure(CalculusStatus::non_finite, CalculusReason::non_finite_result, CalculusPhase::quadrature, "composite Newton-Cotes produced a non-finite result", options, evaluations);
    return detail::integration_success(result, options, evaluations);
}

inline IntegrationResult integrate_romberg(const Function& f, double a, double b, IntegrationOptions options = {}) {
    if (!f) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::invalid_callable, CalculusPhase::setup, "integration function is empty", options);
    if (!options.valid()) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::invalid_options, CalculusPhase::setup, "integration options are invalid", options);
    if (!detail::finite(a) || !detail::finite(b)) return detail::integration_failure(CalculusStatus::invalid_input, CalculusReason::non_finite_input, CalculusPhase::setup, "integration interval must be finite", options);
    if (a == b) return detail::integration_success(0.0, options, 0, 0, 0.0, {{0.0}});

    const double sign = detail::oriented_sign(a, b);
    auto [left, right] = detail::ordered_interval(a, b);
    const double width = right - left;
    std::vector<std::vector<double>> table(static_cast<std::size_t>(std::max(options.max_iterations + 1, 1)),
                                           std::vector<double>(4, 0.0));
    int evaluations = 0;
    CalculusDiagnostic failure;
    const double fa = detail::checked_eval(f, left, failure, evaluations, CalculusPhase::refinement, 0);
    if (failure.status != CalculusStatus::success) return detail::integration_failure(failure.status, failure.reason, failure.phase, failure.message, options, evaluations, 0, detail::nan(), table);
    const double fb = detail::checked_eval(f, right, failure, evaluations, CalculusPhase::refinement, 0);
    if (failure.status != CalculusStatus::success) return detail::integration_failure(failure.status, failure.reason, failure.phase, failure.message, options, evaluations, 0, detail::nan(), table);
    table[0][0] = width * 0.5 * (fa + fb);
    if (options.max_iterations == 0) {
        return detail::integration_failure(CalculusStatus::not_converged, CalculusReason::maximum_iterations, CalculusPhase::convergence, "Romberg reached maximum iterations", options, evaluations, 0, sign * table[0][0], table);
    }

    double best_error = std::numeric_limits<double>::infinity();
    for (int iter = 1; iter <= options.max_iterations; ++iter) {
        const double panels = std::pow(2.0, iter - 1);
        double new_sum = 0.0;
        for (int i = 1; i <= static_cast<int>(panels); ++i) {
            const double x = left + (2.0 * i - 1.0) * width / (2.0 * panels);
            const double y = detail::checked_eval(f, x, failure, evaluations, CalculusPhase::refinement, iter);
            if (failure.status != CalculusStatus::success) return detail::integration_failure(failure.status, failure.reason, failure.phase, failure.message, options, evaluations, iter, detail::nan(), table);
            new_sum += y;
        }
        table[static_cast<std::size_t>(iter)][0] = 0.5 * table[static_cast<std::size_t>(iter - 1)][0] + width * new_sum / (2.0 * panels);
        if (!detail::finite(table[static_cast<std::size_t>(iter)][0])) return detail::integration_failure(CalculusStatus::non_finite, CalculusReason::non_finite_result, CalculusPhase::refinement, "Romberg trapezoid refinement became non-finite", options, evaluations, iter, detail::nan(), table);
        best_error = std::abs(table[static_cast<std::size_t>(iter)][0] - table[static_cast<std::size_t>(iter - 1)][0]);

        const int max_col = std::min(iter, 3);
        for (int col = 1; col <= max_col; ++col) {
            const double factor = std::pow(4.0, col);
            table[static_cast<std::size_t>(iter)][static_cast<std::size_t>(col)] =
                (factor * table[static_cast<std::size_t>(iter)][static_cast<std::size_t>(col - 1)] -
                 table[static_cast<std::size_t>(iter - 1)][static_cast<std::size_t>(col - 1)]) / (factor - 1.0);
            const double horizontal = std::abs(table[static_cast<std::size_t>(iter)][static_cast<std::size_t>(col)] -
                                               table[static_cast<std::size_t>(iter)][static_cast<std::size_t>(col - 1)]);
            best_error = std::min(best_error, horizontal);
            if (iter > col) {
                const double vertical = std::abs(table[static_cast<std::size_t>(iter)][static_cast<std::size_t>(col)] -
                                                 table[static_cast<std::size_t>(iter - 1)][static_cast<std::size_t>(col)]);
                best_error = std::min(best_error, vertical);
            }
            if (!detail::finite(table[static_cast<std::size_t>(iter)][static_cast<std::size_t>(col)])) return detail::integration_failure(CalculusStatus::non_finite, CalculusReason::non_finite_result, CalculusPhase::refinement, "Romberg extrapolation became non-finite", options, evaluations, iter, detail::nan(), table);
            if (best_error <= options.tolerance) {
                return detail::integration_success(sign * table[static_cast<std::size_t>(iter)][static_cast<std::size_t>(col)],
                                                   options, evaluations, iter, best_error, table);
            }
        }
        if (best_error <= options.tolerance) {
            return detail::integration_success(sign * table[static_cast<std::size_t>(iter)][0], options, evaluations, iter, best_error, table);
        }
    }
    const int last = options.max_iterations;
    return detail::integration_failure(CalculusStatus::not_converged, CalculusReason::maximum_iterations, CalculusPhase::convergence,
                                       "Romberg reached maximum iterations", options, evaluations, last,
                                       sign * table[static_cast<std::size_t>(last)][std::min(last, 3)], table);
}

} // namespace MatCal::Calculus

#endif
