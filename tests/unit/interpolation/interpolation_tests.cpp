#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Insert.hpp"
#include "MatCal/Interpolation/CubicSpline.hpp"
#include "MatCal/Interpolation/LinearInterpolator.hpp"
#include "linalg/linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

void linear_tests() {
    using MatCal::Interpolation::ExtrapolationPolicy;
    using MatCal::Interpolation::LinearInterpolator;

    LinearInterpolator line({0.0, 1.0, 3.0}, {1.0, 3.0, 7.0});
    expect_near(line.evaluate(0.0), 1.0, 1.0e-12, "linear hits first node");
    expect_near(line.evaluate(1.0), 3.0, 1.0e-12, "linear hits middle node");
    expect_near(line.evaluate(2.0), 5.0, 1.0e-12, "linear midpoint");
    expect_throw([&] { (void)line.evaluate(-1.0); }, "linear default rejects left extrapolation");

    LinearInterpolator extrap({0.0, 1.0}, {1.0, 3.0}, ExtrapolationPolicy::extrapolate);
    expect_near(extrap.evaluate(-1.0), -1.0, 1.0e-12, "linear explicit extrapolation");

    LinearInterpolator clamp({0.0, 1.0}, {1.0, 3.0}, ExtrapolationPolicy::clamp);
    expect_near(clamp.evaluate(2.0), 3.0, 1.0e-12, "linear explicit clamp");

    expect_throw([] { (void)LinearInterpolator(std::vector<double>{}, std::vector<double>{}); }, "linear rejects empty");
    expect_throw([] { (void)LinearInterpolator(std::vector<double>{1.0}, std::vector<double>{2.0}); }, "linear rejects single node");
    expect_throw([] { (void)LinearInterpolator({0.0, 0.0}, {1.0, 2.0}); }, "linear rejects duplicate x");
    expect_throw([] { (void)LinearInterpolator({1.0, 0.0}, {1.0, 2.0}); }, "linear rejects unsorted x");
    expect_throw([] { (void)LinearInterpolator({0.0, 1.0}, {1.0, std::numeric_limits<double>::infinity()}); },
                 "linear rejects non-finite node");
}

void spline_tests() {
    using MatCal::Interpolation::CubicSpline;
    using MatCal::Interpolation::ExtrapolationPolicy;

    CubicSpline constant({0.0, 1.0, 2.0}, {4.0, 4.0, 4.0});
    expect_near(constant.evaluate(0.5), 4.0, 1.0e-12, "constant spline");
    expect_near(constant.second_derivatives()[1], 0.0, 1.0e-12, "constant second derivative");

    CubicSpline linear({0.0, 2.0}, {1.0, 5.0});
    expect_near(linear.evaluate(1.0), 3.0, 1.0e-12, "two-node spline is linear");
    expect_near(linear.second_derivatives()[0], 0.0, 1.0e-12, "two-node natural left");
    expect_near(linear.second_derivatives()[1], 0.0, 1.0e-12, "two-node natural right");

    CubicSpline known({0.0, 1.0, 2.0}, {0.0, 1.0, 0.0});
    expect_near(known.second_derivatives()[0], 0.0, 1.0e-12, "natural left endpoint");
    expect_near(known.second_derivatives()[1], -3.0, 1.0e-12, "known middle second derivative");
    expect_near(known.second_derivatives()[2], 0.0, 1.0e-12, "natural right endpoint");
    expect_near(known.evaluate(0.5), 0.6875, 1.0e-12, "known natural spline left interval");
    expect_near(known.evaluate(1.5), 0.6875, 1.0e-12, "known natural spline right interval");
    expect_near(known.evaluate(1.0), 1.0, 1.0e-12, "spline hits knot");

    const double eps = 1.0e-6;
    expect_near(known.evaluate(1.0 - eps), known.evaluate(1.0 + eps), 1.0e-5, "spline C0 continuity");
    expect_near(known.derivative(1.0 - eps), known.derivative(1.0 + eps), 1.0e-5, "spline C1 continuity");

    CubicSpline extreme({0.0, 1.0e-6, 1.0}, {0.0, 1.0e-6, 1.0});
    expect_near(extreme.evaluate(0.5), 0.5, 1.0e-10, "spline handles uneven spacing on linear data");

    expect_throw([] { (void)CubicSpline({0.0}, {1.0}); }, "spline rejects single node");
    expect_throw([] { (void)CubicSpline({0.0, 0.0}, {1.0, 2.0}); }, "spline rejects duplicate x");
    expect_throw([] { (void)CubicSpline({0.0, 1.0}, {1.0, std::numeric_limits<double>::quiet_NaN()}); },
                 "spline rejects non-finite node");
    expect_throw([&] { (void)known.evaluate(-0.25); }, "spline default rejects extrapolation");

    CubicSpline extrap({0.0, 1.0}, {0.0, 2.0}, ExtrapolationPolicy::extrapolate);
    expect_near(extrap.evaluate(2.0), 4.0, 1.0e-12, "spline explicit extrapolation");
}

void legacy_differential_tests() {
    using MatCal::Algorithm::Insert::CubicSpline;
    using MatCal::Algorithm::Insert::LinearInsert;
    using MatCal::Interpolation::ExtrapolationPolicy;

    std::vector<double> xs{0.0, 1.0, 3.0};
    std::vector<double> ys{1.0, 3.0, 7.0};
    MatCal::Interpolation::LinearInterpolator core(xs, ys, ExtrapolationPolicy::extrapolate);
    LinearInsert legacy(xs, ys);
    for (double x : {-1.0, 0.0, 0.5, 2.0, 4.0}) {
        expect_near(legacy.calculate(x), core.evaluate(x), 1.0e-12, "legacy LinearInsert delegates core");
    }
    expect_true(legacy.getXs().size() == xs.size(), "legacy LinearInsert preserves getXs");
    expect_true(legacy.getYs().size() == ys.size(), "legacy LinearInsert preserves getYs");

    std::vector<double> sx{0.0, 1.0, 2.0};
    std::vector<double> sy{0.0, 1.0, 0.0};
    MatCal::Interpolation::CubicSpline spline_core(sx, sy, ExtrapolationPolicy::extrapolate);
    CubicSpline spline_legacy(sx, sy);
    for (double x : {-0.5, 0.0, 0.5, 1.0, 1.5, 2.5}) {
        expect_near(spline_legacy.calculate(x), spline_core.evaluate(x), 1.0e-12, "legacy CubicSpline delegates core");
    }
    expect_true(spline_legacy.getXs().size() == sx.size(), "legacy CubicSpline preserves getXs");
    expect_true(spline_legacy.getYs().size() == sy.size(), "legacy CubicSpline preserves getYs");
    expect_near(spline_legacy.getM()[1], -3.0, 1.0e-12, "legacy CubicSpline preserves getM");
}

} // namespace

int main() {
    linear_tests();
    spline_tests();
    legacy_differential_tests();
    return finish("MatCal::Interpolation tests");
}
