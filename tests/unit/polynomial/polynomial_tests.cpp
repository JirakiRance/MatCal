#include <functional>
#include <limits>

#include "MatCal/Polynomial/Polynomial.hpp"
#include "QinJiuShao.hpp"
#include "linalg/linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

void oracle_tests() {
    using MatCal::Polynomial::Polynomial;

    Polynomial p = Polynomial::from_terms({{3, 2.0}, {1, -4.0}, {0, 1.0}});
    expect_true(p.degree() == 3, "polynomial degree");
    expect_near(p.coefficient(3), 2.0, 1e-12, "coefficient degree 3");
    expect_near(p.coefficient(2), 0.0, 1e-12, "missing coefficient is zero");
    expect_near(p.evaluate(2.0), 9.0, 1e-12, "Horner evaluation");

    Polynomial derivative = p.derivative();
    expect_near(derivative.evaluate(2.0), 20.0, 1e-12, "analytic derivative");

    Polynomial integral = p.integral(5.0);
    expect_near(integral.coefficient(0), 5.0, 1e-12, "integral constant");
    expect_near(p.definite_integral(0.0, 1.0), -0.5, 1e-12, "definite integral");

    Polynomial q = Polynomial::from_terms({{1, 1.0}, {0, 2.0}});
    expect_near((p + q).evaluate(2.0), 13.0, 1e-12, "polynomial addition");
    expect_near((p - q).evaluate(2.0), 5.0, 1e-12, "polynomial subtraction");
    expect_near((p * q).evaluate(2.0), 36.0, 1e-12, "polynomial multiplication");
    expect_near((p * 2.0).evaluate(2.0), 18.0, 1e-12, "scalar multiplication");
    expect_near((p / 2.0).evaluate(2.0), 4.5, 1e-12, "scalar division");

    Polynomial zero;
    expect_true(zero.is_zero(), "zero polynomial is explicit");
    expect_true(zero.degree() == 0, "zero polynomial degree contract");
    expect_near(zero.evaluate(12.0), 0.0, 1e-12, "zero polynomial evaluates to zero");
    expect_true((Polynomial{1.0, -1.0} + Polynomial{-1.0, 1.0}).is_zero(), "canonical zero trimming");

    Polynomial sparse = Polynomial::from_terms({{12, 3.0}, {0, -1.0}});
    expect_true(sparse.degree() == 12, "high sparse term degree");
    expect_near(sparse.evaluate(2.0), 12287.0, 1e-12, "high sparse term evaluation");
    expect_throw([] { (void)Polynomial::from_terms({{Polynomial::max_dense_degree + 1, 1.0}}); },
                 "huge sparse degree rejected before dense allocation");

    std::function<double(double)> callable;
    {
        Polynomial local = Polynomial::from_terms({{2, 3.0}, {0, 1.0}});
        callable = local.to_function();
    }
    expect_near(callable(2.0), 13.0, 1e-12, "owning callable outlives polynomial");
}

void negative_tests() {
    using MatCal::Polynomial::Polynomial;

    expect_throw([] { (void)Polynomial::from_terms({{0, std::numeric_limits<double>::infinity()}}); },
                 "non-finite coefficient rejected");
    expect_throw([] { (void)Polynomial{1.0}.evaluate(std::numeric_limits<double>::quiet_NaN()); },
                 "non-finite evaluation point rejected");
    expect_throw([] { (void)(Polynomial{1.0} / 0.0); }, "division by zero rejected");
    expect_throw([] { (void)Polynomial{std::numeric_limits<double>::max(), std::numeric_limits<double>::max()}.evaluate(2.0); },
                 "non-finite evaluation output rejected");
}

void legacy_differential_tests() {
    using MatCal::Polynomial::Polynomial;
    using MatCal::Utils::QinJiuShao;

    Polynomial core = Polynomial::from_terms({{4, 1.5}, {2, -3.0}, {0, 2.0}});
    QinJiuShao legacy({{4, 1.5}, {2, -3.0}, {0, 2.0}});

    for (double x : {-2.0, -0.5, 0.0, 1.25, 3.0}) {
        expect_near(legacy.calculate(x), core.evaluate(x), 1e-12, "legacy calculate delegates core");
        expect_near(legacy.derivative().calculate(x), core.derivative().evaluate(x), 1e-12, "legacy derivative delegates core");
        expect_near(legacy.integral(7.0).calculate(x), core.integral(7.0).evaluate(x), 1e-12, "legacy integral delegates core");
    }

    QinJiuShao other({{1, 2.0}, {0, -1.0}});
    Polynomial other_core = Polynomial::from_terms({{1, 2.0}, {0, -1.0}});
    for (double x : {-1.0, 0.25, 2.0}) {
        expect_near((legacy + other).calculate(x), (core + other_core).evaluate(x), 1e-12, "legacy addition delegates core");
        expect_near((legacy - other).calculate(x), (core - other_core).evaluate(x), 1e-12, "legacy subtraction delegates core");
        expect_near((legacy * other).calculate(x), (core * other_core).evaluate(x), 1e-12, "legacy multiplication delegates core");
        expect_near((legacy * 3.0).calculate(x), (core * 3.0).evaluate(x), 1e-12, "legacy scalar multiply delegates core");
    }

    expect_near(legacy.definiteIntegral(-1.0, 2.0), core.definite_integral(-1.0, 2.0), 1e-12,
                "legacy definite integral delegates core");
    QinJiuShao product_other({{1, -1.0}, {0, 2.0}});
    expect_near(legacy.product(product_other, -1.0, 1.0),
                (core * Polynomial::from_terms({{1, -1.0}, {0, 2.0}})).definite_integral(-1.0, 1.0),
                1e-12,
                "legacy product uses delegated polynomial operations");
    auto f = legacy.toFunction();
    expect_near(f(2.0), core.evaluate(2.0), 1e-12, "legacy owning callable delegates core");
}

} // namespace

int main() {
    oracle_tests();
    negative_tests();
    legacy_differential_tests();
    return finish("MatCal::Polynomial tests");
}
