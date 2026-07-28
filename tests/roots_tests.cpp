#include <cmath>
#include <limits>

#include "Iteration.hpp"
#include "MatCal/Roots/Roots.hpp"
#include "linalg/linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

MatCal::Roots::RootOptions strict_options(int max_iterations = 100) {
    MatCal::Roots::RootOptions options;
    options.absolute_tolerance = 1.0e-12;
    options.relative_tolerance = 1.0e-12;
    options.max_iterations = max_iterations;
    return options;
}

void oracle_tests() {
    using namespace MatCal::Roots;
    auto square_minus_two = [](double x) { return x * x - 2.0; };
    auto two_x = [](double x) { return 2.0 * x; };

    RootResult bisection = solve_bisection(square_minus_two, 0.0, 2.0, strict_options());
    expect_true(bisection.success(), "bisection converges");
    expect_near(bisection.value, std::sqrt(2.0), 1.0e-10, "bisection sqrt2");
    expect_true(bisection.metrics.final_step <= bisection.metrics.tolerance_used, "bisection interval-step contract");

    RootResult newton = solve_newton(square_minus_two, two_x, 1.0, strict_options());
    expect_true(newton.success(), "Newton converges");
    expect_near(newton.value, std::sqrt(2.0), 1.0e-12, "Newton sqrt2");

    RootResult downhill = solve_downhill_newton(square_minus_two, two_x, 10.0, strict_options());
    expect_true(downhill.success(), "downhill Newton converges");
    expect_near(downhill.value, std::sqrt(2.0), 1.0e-10, "downhill Newton sqrt2");

    auto fixed_point = [](double x) { return std::cos(x); };
    RootResult picard = solve_picard(fixed_point, 0.5, strict_options(200));
    expect_true(picard.success(), "Picard converges on cos fixed point");
    expect_near(picard.value, 0.7390851332151607, 1.0e-10, "Picard cos fixed point");

    MatCal::Roots::RootOptions aitken_options = strict_options(200);
    aitken_options.absolute_tolerance = 1.0e-10;
    RootResult aitken = solve_picard_aitken(fixed_point, 0.5, aitken_options);
    expect_true(aitken.success(), "Aitken converges on cos fixed point");
    expect_near(aitken.value, 0.7390851332151607, 1.0e-8, "Aitken cos fixed point");

    RootResult secant_two = solve_secant_two_point(square_minus_two, 1.0, 2.0, strict_options());
    expect_true(secant_two.success(), "two-point secant converges");
    expect_near(secant_two.value, std::sqrt(2.0), 1.0e-10, "two-point secant sqrt2");

    RootResult secant_one = solve_secant_one_point(square_minus_two, 1.5, 1.0e-4, strict_options());
    expect_true(secant_one.success(), "one-point secant converges");
    expect_near(secant_one.value, std::sqrt(2.0), 1.0e-10, "one-point secant sqrt2");

    auto cubic = [](double x) { return x * x * x; };
    auto cubic_derivative = [](double x) { return 3.0 * x * x; };
    RootResult cubic_root = solve_newton(cubic, cubic_derivative, 0.1, strict_options(200));
    expect_true(cubic_root.success(), "Newton handles flat cubic root when not starting at zero derivative");
    expect_near(cubic_root.value, 0.0, 1.0e-4, "flat cubic root");
}

void failure_tests() {
    using namespace MatCal::Roots;
    auto square_plus_one = [](double x) { return x * x + 1.0; };
    auto square = [](double x) { return x * x; };
    auto zero = [](double) { return 0.0; };
    auto non_finite = [](double) { return std::numeric_limits<double>::infinity(); };

    RootResult no_bracket = solve_bisection(square_plus_one, -1.0, 1.0, strict_options());
    expect_true(no_bracket.diagnostic.status == RootStatus::bracket_error, "bisection rejects missing sign change");
    expect_true(no_bracket.diagnostic.reason == RootReason::no_sign_change, "bisection reason no sign change");

    RootResult zero_derivative = solve_newton(square, zero, 1.0, strict_options());
    expect_true(zero_derivative.diagnostic.status == RootStatus::derivative_zero, "Newton reports zero derivative");

    RootResult secant_denominator = solve_secant_two_point(square_plus_one, 1.0, -1.0, strict_options());
    expect_true(secant_denominator.diagnostic.status == RootStatus::denominator_zero, "secant reports zero denominator");

    RootResult non_finite_result = solve_bisection(non_finite, -1.0, 1.0, strict_options());
    expect_true(non_finite_result.diagnostic.status == RootStatus::non_finite, "non-finite callable rejected");

    RootOptions zero_iteration = strict_options(0);
    RootResult max_zero = solve_newton(square_plus_one, [](double x) { return 2.0 * x; }, 1.0, zero_iteration);
    expect_true(!max_zero.success(), "max_iterations=0 does not report success");
    expect_true(max_zero.diagnostic.status == RootStatus::not_converged, "max_iterations=0 reports not converged");

    RootOptions bad = strict_options();
    bad.absolute_tolerance = -1.0;
    RootResult invalid = solve_picard([](double x) { return x; }, 1.0, bad);
    expect_true(invalid.diagnostic.status == RootStatus::invalid_input, "invalid options rejected");

    auto scaled_small = [](double x) { return 1.0e-20 * (x - 3.0); };
    auto scaled_large = [](double x) { return 1.0e20 * (x - 3.0); };
    RootResult small = solve_bisection(scaled_small, 1.0, 5.0, strict_options());
    RootResult large = solve_bisection(scaled_large, 1.0, 5.0, strict_options());
    expect_true(small.success(), "small-scale root succeeds");
    expect_true(large.success(), "large-scale root succeeds");
    expect_near(small.value, 3.0, 1.0e-10, "small-scale root value");
    expect_near(large.value, 3.0, 1.0e-10, "large-scale root value");
}

void legacy_differential_tests() {
    using MatCal::Algorithm::Iteration::Bisection;
    using MatCal::Algorithm::Iteration::Newton;
    using MatCal::Algorithm::Iteration::Picard;
    using MatCal::Algorithm::Iteration::Secant;
    auto square_minus_two = [](double x) { return x * x - 2.0; };
    auto two_x = [](double x) { return 2.0 * x; };

    auto core = MatCal::Roots::solve_bisection(square_minus_two, 0.0, 2.0, strict_options());
    auto legacy = Bisection::solveDetailed(square_minus_two, 0.0, 2.0, 1.0e-12, 100);
    expect_true(legacy.converged == core.converged, "legacy bisection convergence matches");
    expect_near(legacy.root, core.value, 1.0e-12, "legacy bisection delegates core");

    auto ncore = MatCal::Roots::solve_newton(square_minus_two, two_x, 1.0, strict_options());
    auto nlegacy = Newton::solve(square_minus_two, two_x, 1.0, 1.0e-12, 100);
    expect_true(nlegacy.converged == ncore.converged, "legacy Newton convergence matches");
    expect_near(nlegacy.root, ncore.value, 1.0e-12, "legacy Newton delegates core");

    auto pcore = MatCal::Roots::solve_picard([](double x) { return std::cos(x); }, 0.5, strict_options(200));
    auto plegacy = Picard::solveDetailed([](double x) { return std::cos(x); }, 0.5, 1.0e-12, 200);
    expect_true(plegacy.converged == pcore.converged, "legacy Picard convergence matches");
    expect_near(plegacy.root, pcore.value, 1.0e-12, "legacy Picard delegates core");

    auto score = MatCal::Roots::solve_secant_two_point(square_minus_two, 1.0, 2.0, strict_options());
    auto slegacy = Secant::solve_two_point(square_minus_two, 1.0, 2.0, 1.0e-12, 100);
    expect_true(slegacy.converged == score.converged, "legacy secant convergence matches");
    expect_near(slegacy.root, score.value, 1.0e-12, "legacy secant delegates core");
}

} // namespace

int main() {
    oracle_tests();
    failure_tests();
    legacy_differential_tests();
    return finish("MatCal::Roots tests");
}
