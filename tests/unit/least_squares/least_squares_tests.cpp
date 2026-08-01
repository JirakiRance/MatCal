#include <limits>
#include <vector>

#include "test_support.hpp"
#include "Basics.hpp"
#include "MatCal/LeastSquares/LeastSquares.hpp"

namespace {

using matcal_test::expect_near;
using matcal_test::expect_true;
using matcal_test::finish;

void exact_polynomial_fits() {
    std::vector<double> x{-2.0, -1.0, 0.0, 1.0, 2.0};

    std::vector<double> constant{3.0, 3.0, 3.0, 3.0, 3.0};
    auto c = MatCal::LeastSquares::fit_polynomial_degree(0, x, constant);
    expect_true(c.success(), "constant least-squares fit succeeds");
    expect_near(c.polynomial.evaluate(1.5), 3.0, 1.0e-12, "constant fit value");

    std::vector<double> line;
    std::vector<double> quadratic;
    for (double value : x) {
        line.push_back(2.0 * value + 1.0);
        quadratic.push_back(3.0 * value * value - 2.0 * value + 1.0);
    }
    auto l = MatCal::LeastSquares::fit_polynomial_degree(1, x, line);
    auto q = MatCal::LeastSquares::fit_polynomial_degree(2, x, quadratic);
    expect_true(l.success(), "line least-squares fit succeeds");
    expect_true(q.success(), "quadratic least-squares fit succeeds");
    expect_near(l.polynomial.evaluate(1.25), 3.5, 1.0e-11, "line fit value");
    expect_near(q.polynomial.evaluate(1.25), 3.1875, 1.0e-10, "quadratic fit value");
}

void weighted_and_selected_terms() {
    std::vector<double> x{-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<double> y;
    for (double value : x) {
        y.push_back(1.0 + 3.0 * value * value);
    }
    std::vector<double> weights(x.size(), 1.0);
    auto selected = MatCal::LeastSquares::fit_polynomial_selected(2, x, y, weights, {true, false, true});
    expect_true(selected.success(), "selected-term least-squares fit succeeds");
    expect_true(selected.coefficients.size() == 2, "selected-term result only reports selected coefficients");
    expect_near(selected.polynomial.evaluate(1.5), 7.75, 1.0e-10, "selected-term fit value");

    y[2] = 2.0;
    weights[2] = 100.0;
    auto weighted = MatCal::LeastSquares::fit_polynomial_degree(1, x, y, weights);
    expect_true(weighted.success(), "weighted least-squares fit succeeds");
    expect_true(weighted.metrics.residual_norm_inf > 0.0, "weighted noisy fit reports residual");
}

void invalid_and_rank_deficient_inputs() {
    auto rank = MatCal::LeastSquares::fit_polynomial_degree(1, {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0});
    expect_true(rank.status == MatCal::LeastSquares::LeastSquaresStatus::rank_deficient,
                "rank deficient normal equations are detected");

    auto mismatch = MatCal::LeastSquares::fit_polynomial_degree(1, {0.0, 1.0}, {1.0});
    expect_true(mismatch.diagnostic.reason == MatCal::LeastSquares::LeastSquaresReason::size_mismatch,
                "least-squares size mismatch is structured");

    auto bad_weight = MatCal::LeastSquares::fit_polynomial_degree(1, {0.0, 1.0}, {1.0, 2.0}, {1.0, 0.0});
    expect_true(bad_weight.diagnostic.reason == MatCal::LeastSquares::LeastSquaresReason::invalid_weight,
                "least-squares rejects non-positive weights");

    auto nonfinite = MatCal::LeastSquares::fit_polynomial_degree(
        1, {0.0, std::numeric_limits<double>::infinity()}, {1.0, 2.0});
    expect_true(nonfinite.status == MatCal::LeastSquares::LeastSquaresStatus::non_finite_input,
                "least-squares rejects non-finite samples");
}

void legacy_least_square_delegates_core() {
    using Legacy = MatCal::Algorithm::Basics::Least_Square;
    std::vector<double> x{-1.0, 0.0, 1.0, 2.0};
    std::vector<double> y;
    for (double value : x) {
        y.push_back(2.0 * value + 1.0);
    }

    auto legacy = Legacy::solve(1, x, y);
    auto core = MatCal::LeastSquares::fit_polynomial_degree(1, x, y);
    expect_true(legacy.msg == "success!", "legacy least-squares success maps from core");
    expect_true(core.success(), "core least-squares success for differential");
    expect_near(legacy.poly.calculate(1.5), core.polynomial.evaluate(1.5), 1.0e-10,
                "legacy/core least-squares differential");

    std::vector<bool> selects{true, false, true};
    std::vector<double> q;
    for (double value : x) {
        q.push_back(1.0 + 3.0 * value * value);
    }
    auto legacy_selected = Legacy::solve(2, x, q, selects);
    auto core_selected = MatCal::LeastSquares::fit_polynomial_selected(2, x, q, std::vector<double>(x.size(), 1.0), selects);
    expect_near(legacy_selected.poly.calculate(1.5), core_selected.polynomial.evaluate(1.5), 1.0e-10,
                "legacy/core selected-term least-squares differential");
}

} // namespace

int main() {
    exact_polynomial_fits();
    weighted_and_selected_terms();
    invalid_and_rank_deficient_inputs();
    legacy_least_square_delegates_core();
    return finish("MatCal::LeastSquares tests");
}
