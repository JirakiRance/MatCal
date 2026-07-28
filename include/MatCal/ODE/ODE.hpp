#ifndef MATCAL_ODE_ODE_HPP
#define MATCAL_ODE_ODE_HPP

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace MatCal::ODE {

enum class OdeStatus {
    success,
    invalid_input,
    non_finite,
    size_mismatch
};

enum class OdeReason {
    none,
    invalid_rhs,
    invalid_state,
    invalid_time,
    invalid_step,
    invalid_count,
    rhs_size_mismatch,
    non_finite_input,
    non_finite_rhs,
    non_finite_state
};

enum class OdePhase {
    setup,
    rhs,
    step,
    integration
};

struct OdeDiagnostic {
    OdeStatus status = OdeStatus::success;
    OdeReason reason = OdeReason::none;
    OdePhase phase = OdePhase::setup;
    int step_index = 0;
    std::size_t component = 0;
    std::string message;
};

struct OdeOptions {
    bool allow_negative_dt = true;
    int step_count = 1;

    bool valid() const noexcept {
        return step_count >= 0;
    }
};

struct OdeMetrics {
    int steps = 0;
    int rhs_evaluations = 0;
};

struct OdeStepResult {
    std::vector<double> state;
    double t = 0.0;
    bool completed = false;
    OdeDiagnostic diagnostic;
    OdeMetrics metrics;

    bool success() const noexcept {
        return completed && diagnostic.status == OdeStatus::success;
    }
};

struct OdeTrajectoryResult {
    std::vector<std::vector<double>> trajectory;
    bool completed = false;
    OdeDiagnostic diagnostic;
    OdeMetrics metrics;

    bool success() const noexcept {
        return completed && diagnostic.status == OdeStatus::success;
    }
};

using Rhs = std::function<std::vector<double>(double t, const std::vector<double>& y)>;

namespace detail {

inline bool finite(double value) noexcept {
    return std::isfinite(value);
}

inline OdeDiagnostic success_diag() {
    return {OdeStatus::success, OdeReason::none, OdePhase::step, 0, 0, "success"};
}

inline OdeDiagnostic failure_diag(OdeStatus status,
                                  OdeReason reason,
                                  OdePhase phase,
                                  const std::string& message,
                                  int step_index = 0,
                                  std::size_t component = 0) {
    return {status, reason, phase, step_index, component, message};
}

inline bool state_finite(const std::vector<double>& state, std::size_t* bad = nullptr) {
    for (std::size_t i = 0; i < state.size(); ++i) {
        if (!finite(state[i])) {
            if (bad) *bad = i;
            return false;
        }
    }
    return true;
}

inline OdeStepResult step_failure(OdeStatus status,
                                  OdeReason reason,
                                  OdePhase phase,
                                  const std::string& message,
                                  double t,
                                  std::vector<double> state = {},
                                  int rhs_evaluations = 0,
                                  int step_index = 0,
                                  std::size_t component = 0) {
    OdeStepResult result;
    result.state = std::move(state);
    result.t = t;
    result.completed = false;
    result.diagnostic = failure_diag(status, reason, phase, message, step_index, component);
    result.metrics.rhs_evaluations = rhs_evaluations;
    return result;
}

inline OdeTrajectoryResult trajectory_failure(OdeStatus status,
                                              OdeReason reason,
                                              OdePhase phase,
                                              const std::string& message,
                                              std::vector<std::vector<double>> trajectory,
                                              int rhs_evaluations,
                                              int steps,
                                              int step_index = 0,
                                              std::size_t component = 0) {
    OdeTrajectoryResult result;
    result.trajectory = std::move(trajectory);
    result.completed = false;
    result.diagnostic = failure_diag(status, reason, phase, message, step_index, component);
    result.metrics.steps = steps;
    result.metrics.rhs_evaluations = rhs_evaluations;
    return result;
}

inline std::vector<double> add_scaled(double t,
                                      const std::vector<double>& y,
                                      const std::vector<double>& k,
                                      double dt_scale) {
    std::vector<double> out(y.size() + 1, 0.0);
    out[0] = t;
    for (std::size_t i = 0; i < y.size(); ++i) {
        out[i + 1] = y[i] + dt_scale * k[i];
    }
    return out;
}

inline std::vector<double> strip_time(const std::vector<double>& row) {
    if (row.empty()) {
        return {};
    }
    return std::vector<double>(row.begin() + 1, row.end());
}

inline std::vector<double> combine_time_state(double t, const std::vector<double>& state) {
    std::vector<double> row(state.size() + 1, 0.0);
    row[0] = t;
    for (std::size_t i = 0; i < state.size(); ++i) row[i + 1] = state[i];
    return row;
}

inline bool validate_common(const Rhs& rhs,
                            double t,
                            double dt,
                            const std::vector<double>& state,
                            const OdeOptions& options,
                            OdeStepResult& failure) {
    if (!rhs) {
        failure = step_failure(OdeStatus::invalid_input, OdeReason::invalid_rhs, OdePhase::setup, "ODE RHS is empty", t, state);
        return false;
    }
    if (!options.valid()) {
        failure = step_failure(OdeStatus::invalid_input, OdeReason::invalid_count, OdePhase::setup, "ODE options are invalid", t, state);
        return false;
    }
    if (state.empty()) {
        failure = step_failure(OdeStatus::invalid_input, OdeReason::invalid_state, OdePhase::setup, "ODE state dimension is zero", t, state);
        return false;
    }
    std::size_t bad = 0;
    if (!finite(t) || !finite(dt)) {
        failure = step_failure(OdeStatus::invalid_input, OdeReason::non_finite_input, OdePhase::setup, "ODE time and dt must be finite", t, state);
        return false;
    }
    if (!options.allow_negative_dt && dt < 0.0) {
        failure = step_failure(OdeStatus::invalid_input, OdeReason::invalid_step, OdePhase::setup, "negative dt is not allowed by options", t, state);
        return false;
    }
    if (!state_finite(state, &bad)) {
        failure = step_failure(OdeStatus::invalid_input, OdeReason::non_finite_input, OdePhase::setup, "ODE state contains a non-finite value", t, state, 0, 0, bad);
        return false;
    }
    return true;
}

inline std::vector<double> eval_rhs(const Rhs& rhs,
                                    double t,
                                    const std::vector<double>& state,
                                    OdeStepResult& failure,
                                    int& evaluations,
                                    int step_index = 0) {
    std::vector<double> values = rhs(t, state);
    ++evaluations;
    if (values.size() != state.size()) {
        failure = step_failure(OdeStatus::size_mismatch, OdeReason::rhs_size_mismatch, OdePhase::rhs,
                               "ODE RHS returned a vector with the wrong size", t, state, evaluations, step_index);
        return {};
    }
    std::size_t bad = 0;
    if (!state_finite(values, &bad)) {
        failure = step_failure(OdeStatus::non_finite, OdeReason::non_finite_rhs, OdePhase::rhs,
                               "ODE RHS returned a non-finite value", t, state, evaluations, step_index, bad);
        return {};
    }
    return values;
}

} // namespace detail

inline OdeStepResult rk4_step(const Rhs& rhs, double t, const std::vector<double>& state, double dt, OdeOptions options = {}) {
    OdeStepResult failure;
    if (!detail::validate_common(rhs, t, dt, state, options, failure)) return failure;
    if (dt == 0.0) {
        return {state, t, true, detail::success_diag(), {0, 0}};
    }

    int evaluations = 0;
    std::vector<double> k1 = detail::eval_rhs(rhs, t, state, failure, evaluations);
    if (failure.diagnostic.status != OdeStatus::success) return failure;

    std::vector<double> tmp_state(state.size());
    for (std::size_t i = 0; i < state.size(); ++i) tmp_state[i] = state[i] + 0.5 * dt * k1[i];
    std::vector<double> k2 = detail::eval_rhs(rhs, t + 0.5 * dt, tmp_state, failure, evaluations);
    if (failure.diagnostic.status != OdeStatus::success) return failure;

    for (std::size_t i = 0; i < state.size(); ++i) tmp_state[i] = state[i] + 0.5 * dt * k2[i];
    std::vector<double> k3 = detail::eval_rhs(rhs, t + 0.5 * dt, tmp_state, failure, evaluations);
    if (failure.diagnostic.status != OdeStatus::success) return failure;

    for (std::size_t i = 0; i < state.size(); ++i) tmp_state[i] = state[i] + dt * k3[i];
    std::vector<double> k4 = detail::eval_rhs(rhs, t + dt, tmp_state, failure, evaluations);
    if (failure.diagnostic.status != OdeStatus::success) return failure;

    std::vector<double> next(state.size());
    const double c = dt / 6.0;
    for (std::size_t i = 0; i < state.size(); ++i) {
        next[i] = state[i] + c * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    std::size_t bad = 0;
    if (!detail::state_finite(next, &bad) || !detail::finite(t + dt)) {
        return detail::step_failure(OdeStatus::non_finite, OdeReason::non_finite_state, OdePhase::step,
                                    "RK4 produced a non-finite state", t + dt, next, evaluations, 0, bad);
    }
    return {next, t + dt, true, detail::success_diag(), {1, evaluations}};
}

inline OdeStepResult euler_step(const Rhs& rhs, double t, const std::vector<double>& state, double dt, OdeOptions options = {}) {
    OdeStepResult failure;
    if (!detail::validate_common(rhs, t, dt, state, options, failure)) return failure;
    int evaluations = 0;
    std::vector<double> k1 = detail::eval_rhs(rhs, t, state, failure, evaluations);
    if (failure.diagnostic.status != OdeStatus::success) return failure;
    std::vector<double> next(state.size());
    for (std::size_t i = 0; i < state.size(); ++i) next[i] = state[i] + dt * k1[i];
    std::size_t bad = 0;
    if (!detail::state_finite(next, &bad)) {
        return detail::step_failure(OdeStatus::non_finite, OdeReason::non_finite_state, OdePhase::step,
                                    "Euler produced a non-finite state", t + dt, next, evaluations, 0, bad);
    }
    return {next, t + dt, true, detail::success_diag(), {1, evaluations}};
}

inline std::pair<OdeStepResult, std::vector<double>> improved_euler_step(const Rhs& rhs, double t, const std::vector<double>& state, double dt, OdeOptions options = {}) {
    OdeStepResult failure;
    if (!detail::validate_common(rhs, t, dt, state, options, failure)) return {failure, {}};
    int evaluations = 0;
    std::vector<double> k1 = detail::eval_rhs(rhs, t, state, failure, evaluations);
    if (failure.diagnostic.status != OdeStatus::success) return {failure, {}};
    std::vector<double> predicted(state.size());
    for (std::size_t i = 0; i < state.size(); ++i) predicted[i] = state[i] + dt * k1[i];
    std::vector<double> k2 = detail::eval_rhs(rhs, t + dt, predicted, failure, evaluations);
    if (failure.diagnostic.status != OdeStatus::success) return {failure, predicted};
    std::vector<double> next(state.size());
    for (std::size_t i = 0; i < state.size(); ++i) next[i] = state[i] + 0.5 * dt * (k1[i] + k2[i]);
    std::size_t bad = 0;
    if (!detail::state_finite(next, &bad)) {
        return {detail::step_failure(OdeStatus::non_finite, OdeReason::non_finite_state, OdePhase::step,
                                     "improved Euler produced a non-finite state", t + dt, next, evaluations, 0, bad),
                predicted};
    }
    return {{next, t + dt, true, detail::success_diag(), {1, evaluations}}, predicted};
}

inline OdeTrajectoryResult integrate_rk4(const Rhs& rhs, double t0, const std::vector<double>& state0, double dt, int count, OdeOptions options = {}) {
    options.step_count = count;
    if (count < 0) {
        return detail::trajectory_failure(OdeStatus::invalid_input, OdeReason::invalid_count, OdePhase::setup,
                                          "RK4 integration count must be non-negative", {}, 0, 0);
    }
    std::vector<std::vector<double>> trajectory;
    trajectory.reserve(static_cast<std::size_t>(count) + 1);
    trajectory.push_back(detail::combine_time_state(t0, state0));
    double t = t0;
    std::vector<double> state = state0;
    int rhs_evaluations = 0;
    for (int i = 1; i <= count; ++i) {
        OdeStepResult step = rk4_step(rhs, t, state, dt, options);
        rhs_evaluations += step.metrics.rhs_evaluations;
        if (!step.success()) {
            return detail::trajectory_failure(step.diagnostic.status, step.diagnostic.reason, step.diagnostic.phase,
                                              step.diagnostic.message, trajectory, rhs_evaluations, i - 1, i, step.diagnostic.component);
        }
        t = step.t;
        state = step.state;
        trajectory.push_back(detail::combine_time_state(t, state));
    }
    return {trajectory, true, detail::success_diag(), {count, rhs_evaluations}};
}

inline OdeTrajectoryResult integrate_euler(const Rhs& rhs, double t0, const std::vector<double>& state0, double dt, int count, OdeOptions options = {}) {
    if (count < 0) {
        return detail::trajectory_failure(OdeStatus::invalid_input, OdeReason::invalid_count, OdePhase::setup,
                                          "Euler integration count must be non-negative", {}, 0, 0);
    }
    std::vector<std::vector<double>> trajectory;
    trajectory.reserve(static_cast<std::size_t>(count) + 1);
    trajectory.push_back(detail::combine_time_state(t0, state0));
    double t = t0;
    std::vector<double> state = state0;
    int rhs_evaluations = 0;
    for (int i = 1; i <= count; ++i) {
        OdeStepResult step = euler_step(rhs, t, state, dt, options);
        rhs_evaluations += step.metrics.rhs_evaluations;
        if (!step.success()) {
            return detail::trajectory_failure(step.diagnostic.status, step.diagnostic.reason, step.diagnostic.phase,
                                              step.diagnostic.message, trajectory, rhs_evaluations, i - 1, i, step.diagnostic.component);
        }
        t = step.t;
        state = step.state;
        trajectory.push_back(detail::combine_time_state(t, state));
    }
    return {trajectory, true, detail::success_diag(), {count, rhs_evaluations}};
}

} // namespace MatCal::ODE

#endif
