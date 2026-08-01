#include <limits>
#include <vector>

#include "test_support.hpp"
#include "Insert.hpp"
#include "MatCal/Interpolation/PolynomialInterpolation.hpp"

namespace {

using matcal_test::expect_near;
using matcal_test::expect_throw;
using matcal_test::expect_true;
using matcal_test::finish;

void lagrange_and_newton_recover_polynomials() {
    std::vector<std::pair<double, double>> data{{-1.0, 2.0}, {0.0, 1.0}, {1.0, 2.0}, {2.0, 5.0}};
    auto lagrange = MatCal::Interpolation::interpolate_lagrange(data);
    auto newton = MatCal::Interpolation::interpolate_newton_divided(data);

    for (const auto& [x, y] : data) {
        expect_near(lagrange.evaluate(x), y, 1.0e-12, "Lagrange hits node");
        expect_near(newton.polynomial.evaluate(x), y, 1.0e-12, "Newton divided hits node");
    }
    expect_near(lagrange.evaluate(1.5), 3.25, 1.0e-12, "Lagrange recovers x^2+1 midpoint");
    expect_near(newton.polynomial.evaluate(1.5), lagrange.evaluate(1.5), 1.0e-12,
                "Newton divided matches Lagrange");

    expect_throw([] {
        (void)MatCal::Interpolation::interpolate_lagrange({{0.0, 1.0}, {0.0, 2.0}});
    }, "Lagrange rejects duplicate x");
    expect_throw([] {
        (void)MatCal::Interpolation::interpolate_newton_divided(
            {{0.0, 1.0}, {std::numeric_limits<double>::quiet_NaN(), 2.0}});
    }, "Newton divided rejects non-finite x");
}

void finite_difference_interpolation_contract() {
    auto finite = MatCal::Interpolation::interpolate_newton_finite(1.0, 0.0, {0.0, 1.0, 4.0, 9.0});
    expect_near(finite.polynomial.evaluate(1.5), 2.25, 1.0e-12, "Newton finite recovers x^2");
    expect_true(finite.table.size() == 4, "Newton finite exposes difference table");
    expect_near(finite.table[2][2], 2.0, 1.0e-12, "Newton finite second difference");
    expect_throw([] {
        (void)MatCal::Interpolation::interpolate_newton_finite(0.0, 0.0, {0.0, 1.0});
    }, "Newton finite rejects zero spacing");
    expect_throw([] {
        (void)MatCal::Interpolation::interpolate_newton_finite(1.0, 0.0, {0.0});
    }, "Newton finite rejects single value");
}

void hermite_matches_values_and_derivatives() {
    std::vector<double> xs{0.0, 1.0};
    std::vector<double> ys{0.0, 1.0};
    std::vector<double> dys{0.0, 2.0};
    auto hermite = MatCal::Interpolation::interpolate_hermite(xs, ys, dys);
    auto derivative = hermite.derivative();
    expect_near(hermite.evaluate(0.0), 0.0, 1.0e-12, "Hermite hits first value");
    expect_near(hermite.evaluate(1.0), 1.0, 1.0e-12, "Hermite hits second value");
    expect_near(derivative.evaluate(0.0), 0.0, 1.0e-12, "Hermite hits first derivative");
    expect_near(derivative.evaluate(1.0), 2.0, 1.0e-12, "Hermite hits second derivative");
    expect_near(hermite.evaluate(0.5), 0.25, 1.0e-12, "Hermite recovers x^2 for two-node data");

    expect_throw([] {
        (void)MatCal::Interpolation::interpolate_hermite({0.0, 0.0}, {0.0, 1.0}, {0.0, 1.0});
    }, "Hermite rejects duplicate x");
}

void legacy_interpolation_delegates_core() {
    using MatCal::Algorithm::Insert::Hermite;
    using MatCal::Algorithm::Insert::LagrangeInsert;
    using MatCal::Algorithm::Insert::NewtonInsert_Finite;
    using MatCal::Algorithm::Insert::NewtonInsert_Quotient;

    std::vector<std::pair<double, double>> data{{0.0, 1.0}, {1.0, 2.0}, {2.0, 5.0}};
    LagrangeInsert legacy_lagrange(data);
    NewtonInsert_Quotient legacy_quotient(data);
    auto core_lagrange = MatCal::Interpolation::interpolate_lagrange(data);
    auto core_quotient = MatCal::Interpolation::interpolate_newton_divided(data);
    expect_near(legacy_lagrange.calculate(1.5), core_lagrange.evaluate(1.5), 1.0e-12,
                "legacy/core Lagrange differential");
    expect_near(legacy_quotient.calculate(1.5), core_quotient.polynomial.evaluate(1.5), 1.0e-12,
                "legacy/core Newton divided differential");
    expect_true(legacy_quotient.getSheet().getRows() == 3, "legacy Newton divided keeps sheet");

    std::vector<double> y{0.0, 1.0, 4.0};
    NewtonInsert_Finite legacy_finite(1.0, 0.0, y);
    auto core_finite = MatCal::Interpolation::interpolate_newton_finite(1.0, 0.0, y);
    expect_near(legacy_finite.calculate(1.5), core_finite.polynomial.evaluate(1.5), 1.0e-12,
                "legacy/core Newton finite differential");
    expect_near(legacy_finite.getSheet().get(2, 2), 2.0, 1.0e-12, "legacy Newton finite keeps difference sheet");

    std::vector<double> xs{0.0, 1.0};
    std::vector<double> ys{0.0, 1.0};
    std::vector<double> dys{0.0, 2.0};
    Hermite legacy_hermite(xs, ys, dys);
    auto core_hermite = MatCal::Interpolation::interpolate_hermite(xs, ys, dys);
    expect_near(legacy_hermite.calculate(0.5), core_hermite.evaluate(0.5), 1.0e-12,
                "legacy/core Hermite differential");
}

} // namespace

int main() {
    lagrange_and_newton_recover_polynomials();
    finite_difference_interpolation_contract();
    hermite_matches_values_and_derivatives();
    legacy_interpolation_delegates_core();
    return finish("remaining interpolation tests");
}
