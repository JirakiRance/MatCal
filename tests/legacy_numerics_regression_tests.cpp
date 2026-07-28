#include <functional>
#include <limits>

#include "Basics.hpp"
#include "QinJiuShao.hpp"
#include "test_support.hpp"

using namespace matcal_test;

namespace {

void instant_regressions() {
    using MatCal::Algorithm::Basics::NumericalIntegration;

    std::function<double(double)> empty;
    expect_throw([&] { (void)NumericalIntegration::Instant(empty, 0.0, 1.0, 0.1); }, "Instant rejects empty callable");
    expect_throw([] { (void)NumericalIntegration::Instant([](double x) { return x; }, 0.0, 1.0, 0.0); }, "Instant rejects zero eps");
    double reverse_constant = NumericalIntegration::Instant([](double) { return 2.0; }, 1.0, 0.0, 0.1);
    expect_near(reverse_constant, -2.0, 1e-12, "Instant preserves signed reversed interval");
    expect_throw([] { (void)NumericalIntegration::Instant([](double x) { return x; }, 0.0, 1.0, std::numeric_limits<double>::infinity()); }, "Instant rejects infinite eps");
    expect_throw([] { (void)NumericalIntegration::Instant([](double x) { return x; }, std::numeric_limits<double>::quiet_NaN(), 1.0, 0.1); }, "Instant rejects NaN interval");

    double constant = NumericalIntegration::Instant([](double) { return 2.0; }, 0.0, 1.0, 0.1);
    expect_near(constant, 2.0, 1e-12, "Instant integrates constant function");

    double linear = NumericalIntegration::Instant([](double x) { return x; }, 0.0, 1.0, 0.1);
    expect_near(linear, 0.45, 1e-12, "Instant preserves legacy left-rectangle rule");
}

void legendre_regressions() {
    using MatCal::Algorithm::Basics::OrthogonalPolynomials;

    auto p0 = OrthogonalPolynomials::Legendre(0);
    auto p1 = OrthogonalPolynomials::Legendre(1);
    auto p2 = OrthogonalPolynomials::Legendre(2);
    auto p3 = OrthogonalPolynomials::Legendre(3);

    for (double x : {-1.0, -0.25, 0.0, 0.5, 1.0}) {
        expect_near(p0.calculate(x), 1.0, 1e-12, "Legendre P0");
        expect_near(p1.calculate(x), x, 1e-12, "Legendre P1");
        expect_near(p2.calculate(x), (3.0 * x * x - 1.0) / 2.0, 1e-12, "Legendre P2");
        expect_near(p3.calculate(x), (5.0 * x * x * x - 3.0 * x) / 2.0, 1e-12, "Legendre P3");
    }
}

void determinant_regressions() {
    using MatCal::Algorithm::Matrix::determinant;
    using MatCal::Utils::Matrix;

    Matrix odd_swap({{0.0, 1.0}, {1.0, 0.0}});
    expect_near(determinant(odd_swap), -1.0, 1e-12, "determinant accounts for odd row swap");
    expect_near(odd_swap.get(0, 0), 0.0, 1e-12, "determinant does not modify input 00");
    expect_near(odd_swap.get(0, 1), 1.0, 1e-12, "determinant does not modify input 01");

    Matrix two_swaps({{0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}});
    expect_near(determinant(two_swaps), 1.0, 1e-12, "determinant accounts for even row swaps");

    Matrix singular({{1.0, 2.0}, {2.0, 4.0}});
    expect_throw([&] { (void)determinant(singular); }, "determinant singular contract throws");
}

void qinjiushao_function_owns_polynomial_state() {
    using MatCal::Utils::QinJiuShao;

    std::function<double(double)> f;
    {
        QinJiuShao p({{2, 3.0}, {1, -2.0}, {0, 1.0}});
        f = p.toFunction();
    }

    expect_near(f(2.0), 9.0, 1e-12, "QinJiuShao toFunction owns copied state");
}

} // namespace

int main() {
    instant_regressions();
    legendre_regressions();
    determinant_regressions();
    qinjiushao_function_owns_polynomial_state();
    return finish("legacy numerics regression tests");
}
