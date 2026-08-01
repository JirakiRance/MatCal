#include <cmath>
#include <limits>

#include "Basics.hpp"
#include "MatCal/Calculus/Calculus.hpp"
#include "linalg/linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

constexpr double pi = 3.14159265358979323846;

void derivative_tests() {
    using namespace MatCal::Calculus;

    auto constant = [](double) { return 5.0; };
    auto square = [](double x) { return x * x; };
    auto sine = [](double x) { return std::sin(x); };
    auto exponential = [](double x) { return std::exp(x); };

    auto c = forward_difference(constant, 2.0, 1.0e-5);
    expect_true(c.success(), "forward derivative constant succeeds");
    expect_near(c.value, 0.0, 1.0e-10, "forward derivative constant");
    expect_true(c.function_evaluations == 2, "forward derivative evaluation count");

    auto sq_forward = forward_difference(square, 3.0, 1.0e-5);
    auto sq_center = central_difference(square, 3.0, 1.0e-5);
    expect_near(sq_forward.value, 6.0, 1.0e-4, "forward derivative x^2");
    expect_near(sq_center.value, 6.0, 1.0e-9, "central derivative x^2");
    expect_true(std::abs(sq_center.value - 6.0) < std::abs(sq_forward.value - 6.0), "central difference improves x^2");

    expect_near(central_difference(sine, 0.3, 1.0e-5).value, std::cos(0.3), 1.0e-10, "central derivative sin");
    expect_near(central_difference(exponential, 0.2, 1.0e-5).value, std::exp(0.2), 1.0e-10, "central derivative exp");

    auto bad_h = forward_difference(square, 1.0, 0.0);
    expect_true(!bad_h.success(), "derivative rejects zero h");
    auto bad_func = central_difference([](double) { return std::numeric_limits<double>::infinity(); }, 1.0, 1.0e-4);
    expect_true(bad_func.diagnostic.status == CalculusStatus::non_finite, "derivative rejects non-finite callable");

    using MatCal::Algorithm::Basics::Derivative;
    expect_near(Derivative::dy_dx(square, 3.0, 1.0e-5), sq_forward.value, 1.0e-12, "legacy dy_dx delegates core");
    expect_near(Derivative::dy_dx_center(square, 3.0, 1.0e-5), sq_center.value, 1.0e-12, "legacy dy_dx_center delegates core");
}

void integration_tests() {
    using namespace MatCal::Calculus;

    auto one = [](double) { return 1.0; };
    auto linear = [](double x) { return x; };
    auto square = [](double x) { return x * x; };
    auto sine = [](double x) { return std::sin(x); };
    auto exponential = [](double x) { return std::exp(x); };

    expect_near(integrate_newton_cotes(one, 2.0, 5.0, 1).value, 3.0, 1.0e-12, "integral constant");
    expect_near(integrate_newton_cotes(linear, 0.0, 1.0, 1).value, 0.5, 1.0e-12, "integral x");
    expect_near(integrate_newton_cotes(square, 0.0, 1.0, 2).value, 1.0 / 3.0, 1.0e-12, "integral x^2");
    expect_near(integrate_composite_newton_cotes(sine, 0.0, pi, 64, 4).value, 2.0, 1.0e-10, "composite integral sin");
    expect_near(integrate_composite_newton_cotes(exponential, 0.0, 1.0, 32, 4).value, std::exp(1.0) - 1.0, 1.0e-10, "composite integral exp");

    expect_near(integrate_newton_cotes(linear, 1.0, 0.0, 1).value, -0.5, 1.0e-12, "reverse interval sign");
    expect_near(integrate_newton_cotes(linear, 2.0, 2.0, 1).value, 0.0, 1.0e-12, "zero interval");
    expect_near(integrate_instant(one, 0.0, 1.0, 0.1).value, 1.0, 1.0e-12, "instant left rectangle constant");

    IntegrationOptions romberg_options;
    romberg_options.tolerance = 1.0e-12;
    romberg_options.max_iterations = 16;
    auto romberg = integrate_romberg(sine, 0.0, pi, romberg_options);
    expect_true(romberg.success(), "Romberg converges for sin");
    expect_near(romberg.value, 2.0, 1.0e-10, "Romberg integral sin");
    expect_true(romberg.metrics.function_evaluations > 2, "Romberg records evaluations");

    IntegrationOptions no_iter;
    no_iter.max_iterations = 0;
    auto not_converged = integrate_romberg(sine, 0.0, pi, no_iter);
    expect_true(!not_converged.success(), "Romberg max_iterations=0 is not success");

    auto bad_func = integrate_newton_cotes([](double) { return std::numeric_limits<double>::quiet_NaN(); }, 0.0, 1.0, 2);
    expect_true(bad_func.diagnostic.status == CalculusStatus::non_finite, "integration rejects non-finite callable");
    auto bad_order = integrate_newton_cotes(one, 0.0, 1.0, 8);
    expect_true(bad_order.diagnostic.reason == CalculusReason::invalid_order, "Newton-Cotes rejects invalid order");
    auto bad_segments = integrate_composite_newton_cotes(one, 0.0, 1.0, 0, 4);
    expect_true(bad_segments.diagnostic.reason == CalculusReason::invalid_segments, "Composite Newton-Cotes rejects invalid segments");

    using MatCal::Algorithm::Basics::NumericalIntegration;
    expect_near(NumericalIntegration::NewtonCotes(square, 0.0, 1.0, 2),
                integrate_newton_cotes(square, 0.0, 1.0, 2).value,
                1.0e-12,
                "legacy NewtonCotes delegates core");
    expect_near(NumericalIntegration::CompositeNewtonCotes(sine, 0.0, pi, 64, 4),
                integrate_composite_newton_cotes(sine, 0.0, pi, 64, 4).value,
                1.0e-12,
                "legacy CompositeNewtonCotes delegates core");
    auto legacy_romberg = NumericalIntegration::Romberg(sine, 0.0, pi, 1.0e-12, 16);
    expect_near(legacy_romberg.first, romberg.value, 1.0e-12, "legacy Romberg delegates core");
    expect_near(NumericalIntegration::Instant(linear, 1.0, 0.0, 0.1), integrate_instant(linear, 1.0, 0.0, 0.1).value, 1.0e-12,
                "legacy Instant delegates reverse interval core");
}

} // namespace

int main() {
    derivative_tests();
    integration_tests();
    return finish("MatCal::Calculus tests");
}
