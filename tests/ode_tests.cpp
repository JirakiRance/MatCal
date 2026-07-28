#include <cmath>
#include <limits>
#include <vector>

#include "Basics.hpp"
#include "MatCal/ODE/ODE.hpp"
#include "linalg/linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

constexpr double pi = 3.14159265358979323846;

void rk4_oracle_tests() {
    using namespace MatCal::ODE;

    auto exp_rhs = [](double, const std::vector<double>& y) {
        return std::vector<double>{y[0]};
    };
    auto decay_rhs = [](double, const std::vector<double>& y) {
        return std::vector<double>{-y[0]};
    };
    auto constant_rhs = [](double, const std::vector<double>&) {
        return std::vector<double>{2.0};
    };

    auto one_step = rk4_step(exp_rhs, 0.0, {1.0}, 0.1);
    expect_true(one_step.success(), "RK4 one step succeeds");
    expect_near(one_step.state[0], std::exp(0.1), 1.0e-7, "RK4 one step y'=y");
    expect_true(one_step.metrics.rhs_evaluations == 4, "RK4 evaluation count");

    auto decay = integrate_rk4(decay_rhs, 0.0, {1.0}, 0.01, 100);
    expect_true(decay.success(), "RK4 decay succeeds");
    expect_near(decay.trajectory.back()[1], std::exp(-1.0), 1.0e-10, "RK4 y'=-y");

    auto constant = integrate_rk4(constant_rhs, 0.0, {3.0}, 0.1, 10);
    expect_near(constant.trajectory.back()[1], 5.0, 1.0e-12, "RK4 constant derivative");

    auto oscillator = [](double, const std::vector<double>& y) {
        return std::vector<double>{y[1], -y[0]};
    };
    auto harmonic = integrate_rk4(oscillator, 0.0, {1.0, 0.0}, 0.01, 628);
    expect_true(harmonic.success(), "RK4 harmonic oscillator succeeds");
    expect_near(harmonic.trajectory.back()[1], std::cos(6.28), 1.0e-7, "RK4 harmonic position");
    expect_near(harmonic.trajectory.back()[2], -std::sin(6.28), 1.0e-7, "RK4 harmonic velocity");

    auto full = integrate_rk4(exp_rhs, 0.0, {1.0}, 0.1, 10);
    auto half = integrate_rk4(exp_rhs, 0.0, {1.0}, 0.05, 20);
    double full_error = std::abs(full.trajectory.back()[1] - std::exp(1.0));
    double half_error = std::abs(half.trajectory.back()[1] - std::exp(1.0));
    expect_true(half_error < full_error / 10.0, "RK4 fourth-order error trend");

    auto backward = integrate_rk4(exp_rhs, 1.0, {std::exp(1.0)}, -0.1, 10);
    expect_true(backward.success(), "RK4 backward integration allowed");
    expect_near(backward.trajectory.back()[1], 1.0, 2.0e-6, "RK4 backward integration");
}

void failure_tests() {
    using namespace MatCal::ODE;
    auto mismatch = [](double, const std::vector<double>&) {
        return std::vector<double>{1.0, 2.0};
    };
    auto non_finite = [](double, const std::vector<double>&) {
        return std::vector<double>{std::numeric_limits<double>::infinity()};
    };

    expect_true(rk4_step({}, 0.0, {1.0}, 0.1).diagnostic.status == OdeStatus::invalid_input, "RK4 rejects empty RHS");
    expect_true(rk4_step(mismatch, 0.0, {1.0}, 0.1).diagnostic.status == OdeStatus::size_mismatch, "RK4 rejects RHS size mismatch");
    expect_true(rk4_step(non_finite, 0.0, {1.0}, 0.1).diagnostic.status == OdeStatus::non_finite, "RK4 rejects non-finite RHS");
    expect_true(rk4_step([](double, const std::vector<double>& y) { return y; }, 0.0, {}, 0.1).diagnostic.status == OdeStatus::invalid_input,
                "RK4 rejects empty state");
}

void legacy_and_pt_differential_tests() {
    using MatCal::Algorithm::Basics::ODE;
    using MatCal::Algorithm::Basics::Integrate::RK4;

    std::vector<std::function<double(std::vector<double>&)>> funcs;
    funcs.push_back([](std::vector<double>& row) { return row[1]; });
    std::vector<double> inits{0.0, 1.0};
    auto legacy = ODE::RungeKutta_44(1, funcs, inits, 0.1, 10);
    auto core = MatCal::ODE::integrate_rk4(
        [](double, const std::vector<double>& y) { return std::vector<double>{y[0]}; },
        0.0,
        {1.0},
        0.1,
        10);
    expect_true(core.success(), "core RK4 trajectory succeeds");
    expect_near(legacy.get(10, 1), core.trajectory.back()[1], 1.0e-12, "legacy RungeKutta_44 delegates core");

    auto simple = ODE::SimpleEuler(1, funcs, inits, 0.1, 2);
    auto euler = MatCal::ODE::integrate_euler(
        [](double, const std::vector<double>& y) { return std::vector<double>{y[0]}; },
        0.0,
        {1.0},
        0.1,
        2);
    expect_near(simple.get(2, 1), euler.trajectory.back()[1], 1.0e-12, "legacy SimpleEuler delegates core");

    std::vector<double> y{1.0};
    std::vector<double> y_out;
    RK4::step([](const std::vector<double>& state, std::vector<double>& dydt) {
        dydt.resize(1);
        dydt[0] = state[0];
    }, y, 0.1, y_out);
    expect_near(y_out[0], std::exp(0.1), 1.0e-7, "PT RK4::step delegates core");

    double theta_out = 0.0;
    double omega_out = 0.0;
    RK4::step2([](double theta, double omega, double& dtheta, double& domega) {
        dtheta = omega;
        domega = -theta;
    }, 1.0, 0.0, 0.1, theta_out, omega_out);
    auto step2_core = MatCal::ODE::rk4_step(
        [](double, const std::vector<double>& state) {
            return std::vector<double>{state[1], -state[0]};
        },
        0.0,
        {1.0, 0.0},
        0.1);
    expect_near(theta_out, step2_core.state[0], 1.0e-12, "PT RK4::step2 theta delegates core");
    expect_near(omega_out, step2_core.state[1], 1.0e-12, "PT RK4::step2 omega delegates core");
}

} // namespace

int main() {
    rk4_oracle_tests();
    failure_tests();
    legacy_and_pt_differential_tests();
    return finish("MatCal::ODE tests");
}
