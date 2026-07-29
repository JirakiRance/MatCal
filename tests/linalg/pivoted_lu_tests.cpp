#include <cmath>
#include <limits>
#include <random>

#include "MatCal/Linalg/DenseSolver.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

MatCal::Linalg::Vector column(const MatCal::Linalg::DenseMatrix& matrix, std::size_t col) {
    MatCal::Linalg::Vector result(matrix.rows());
    for (std::size_t r = 0; r < matrix.rows(); ++r) {
        result[r] = matrix(r, col);
    }
    return result;
}

MatCal::Linalg::DenseMatrix permutation_times_original(
    const MatCal::Linalg::PivotedLuFactorization& factorization) {
    const auto& a = factorization.original_matrix();
    MatCal::Linalg::DenseMatrix result(a.rows(), a.cols());
    const auto& permutation = factorization.row_permutation();
    for (std::size_t r = 0; r < a.rows(); ++r) {
        for (std::size_t c = 0; c < a.cols(); ++c) {
            result(r, c) = a(permutation[r], c);
        }
    }
    return result;
}

MatCal::Linalg::DenseMatrix l_times_u(const MatCal::Linalg::PivotedLuFactorization& factorization) {
    const auto& lu = factorization.lu();
    const auto n = lu.rows();
    MatCal::Linalg::DenseMatrix result(n, n);
    for (std::size_t r = 0; r < n; ++r) {
        for (std::size_t c = 0; c < n; ++c) {
            double sum = 0.0;
            for (std::size_t k = 0; k < n; ++k) {
                const double l = r == k ? 1.0 : (r > k ? lu(r, k) : 0.0);
                const double u = k <= c ? lu(k, c) : 0.0;
                sum += l * u;
            }
            result(r, c) = sum;
        }
    }
    return result;
}

void expect_matrix_near(const MatCal::Linalg::DenseMatrix& actual,
                        const MatCal::Linalg::DenseMatrix& expected,
                        double tolerance,
                        const std::string& label) {
    expect_true(actual.rows() == expected.rows() && actual.cols() == expected.cols(), label + " dimensions");
    for (std::size_t r = 0; r < actual.rows(); ++r) {
        for (std::size_t c = 0; c < actual.cols(); ++c) {
            expect_near(actual(r, c), expected(r, c), tolerance, label + " value");
        }
    }
}

MatCal::Linalg::DenseMatrix scaled_matrix(double scale) {
    return MatCal::Linalg::DenseMatrix{{scale * 4.0, scale * 1.0},
                                       {scale * 2.0, scale * 3.0}};
}

MatCal::Linalg::Vector scaled_rhs(double scale) {
    return MatCal::Linalg::Vector{scale * 1.0, scale * 2.0};
}

void factorization_contracts() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::SolverStatus;
    using MatCal::Linalg::Vector;
    using MatCal::Linalg::factorize_dense_partial_pivot;

    DenseMatrix a{{0.0, 1.0, 1.0}, {2.0, 3.0, 4.0}, {1.0, 1.0, 0.0}};
    const DenseMatrix original = a;
    auto factorization = factorize_dense_partial_pivot(a);
    expect_true(factorization.success(), "pivoted LU factorization succeeds");
    expect_true(factorization.factorization.valid(), "factorization object marked valid");
    expect_true(factorization.factorization.permutation_sign() == -1, "single row swap sign");
    expect_near(a(0, 0), original(0, 0), 0.0, "factorization does not modify input");
    expect_matrix_near(permutation_times_original(factorization.factorization),
                       l_times_u(factorization.factorization),
                       1.0e-12,
                       "PA equals LU");

    Vector rhs{3.0, 20.0, -1.0};
    auto solved = factorization.factorization.solve(rhs);
    expect_true(solved.success(), "factorization solves vector RHS");
    expect_true(solved.metrics.factorization_operation_count == factorization.metrics.factorization_operation_count,
                "solve reports reused factorization count");
    expect_true(solved.metrics.solve_operation_count > 0, "solve operation count recorded");

    auto invalid_rhs = factorization.factorization.solve(DenseMatrix{{1.0}, {std::numeric_limits<double>::quiet_NaN()}, {2.0}});
    expect_true(invalid_rhs.status == SolverStatus::non_finite_input, "multi RHS finite precheck");
    auto solved_after_failure = factorization.factorization.solve(rhs);
    expect_true(solved_after_failure.success(), "failed solve does not invalidate factorization");

    DenseMatrix rhs_multi{{3.0, 6.0}, {20.0, 40.0}, {-1.0, -2.0}};
    auto multi = factorization.factorization.solve(rhs_multi);
    expect_true(multi.success(), "multi RHS solve succeeds");
    expect_near(multi.solution(0, 1), 2.0 * solved.solution[0], 1.0e-12, "multi RHS col1 x0");
    expect_true(multi.metrics.factorization_operation_count == factorization.metrics.factorization_operation_count,
                "multi RHS does not refactor");

    auto one_col0 = MatCal::Linalg::solve_dense_partial_pivot(a, column(rhs_multi, 0));
    auto one_col1 = MatCal::Linalg::solve_dense_partial_pivot(a, column(rhs_multi, 1));
    expect_true(one_col0.success() && one_col1.success(), "one-shot references succeed");
    expect_true(multi.metrics.factorization_operation_count < one_col0.metrics.factorization_operation_count +
                                                        one_col1.metrics.factorization_operation_count,
                "multi RHS factorization count is lower than repeated one-shot");
    expect_near(multi.solution(2, 0), one_col0.solution[2], 1.0e-12, "multi RHS col0 matches one-shot");
    expect_near(multi.solution(2, 1), one_col1.solution[2], 1.0e-12, "multi RHS col1 matches one-shot");
}

void determinant_contracts() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::factorize_dense_partial_pivot;

    auto odd = factorize_dense_partial_pivot(DenseMatrix{{0.0, 1.0}, {1.0, 0.0}});
    expect_true(odd.success(), "odd swap determinant factorization");
    expect_true(odd.factorization.permutation_sign() == -1, "odd swap sign");
    expect_near(odd.factorization.determinant(), -1.0, 0.0, "odd swap determinant");

    auto even = factorize_dense_partial_pivot(DenseMatrix{{0.0, 1.0, 0.0},
                                                         {0.0, 0.0, 1.0},
                                                         {1.0, 0.0, 0.0}});
    expect_true(even.success(), "even swap determinant factorization");
    expect_true(even.factorization.permutation_sign() == 1, "even swap sign");
    expect_near(even.factorization.determinant(), 1.0, 0.0, "even swap determinant");

    auto one = factorize_dense_partial_pivot(DenseMatrix{{5.0}});
    expect_true(one.success(), "1x1 factorization");
    expect_near(one.factorization.determinant(), 5.0, 0.0, "1x1 determinant");

    auto empty = factorize_dense_partial_pivot(DenseMatrix());
    expect_true(empty.success(), "0x0 factorization");
    expect_near(empty.factorization.determinant(), 1.0, 0.0, "0x0 determinant identity");
    auto empty_solve = empty.factorization.solve(DenseMatrix(0, 2));
    expect_true(empty_solve.success() && empty_solve.solution.rows() == 0 && empty_solve.solution.cols() == 2,
                "0x0 factorization solves empty multi RHS");
}

void failure_and_scale_contracts() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::SolverOptions;
    using MatCal::Linalg::SolverStatus;
    using MatCal::Linalg::Vector;
    using MatCal::Linalg::factorize_dense_partial_pivot;
    using MatCal::Linalg::solve_dense_partial_pivot;

    auto singular = factorize_dense_partial_pivot(DenseMatrix{{1.0, 2.0}, {2.0, 4.0}});
    expect_true(singular.status == SolverStatus::singular, "singular factorization rejected");
    expect_true(!singular.factorization.valid(), "singular result has no valid partial factor");

    SolverOptions strict;
    strict.absolute_tolerance = 1.0e-12;
    strict.relative_tolerance = 0.0;
    auto near = factorize_dense_partial_pivot(DenseMatrix{{1.0e-14, 0.0}, {0.0, 1.0}}, strict);
    expect_true(near.status == SolverStatus::singular, "near singular obeys explicit tolerance");

    auto nonfinite = factorize_dense_partial_pivot(
        DenseMatrix{{1.0, std::numeric_limits<double>::infinity()}, {0.0, 1.0}});
    expect_true(nonfinite.status == SolverStatus::non_finite_input, "non-finite matrix rejected");

    for (double scale : {1.0e-20, 1.0, 1.0e20}) {
        auto result = solve_dense_partial_pivot(scaled_matrix(scale), scaled_rhs(scale));
        expect_true(result.success(), "scaled solve success");
        expect_near(result.solution[0], 0.1, 1.0e-12, "scaled solve x0");
        expect_near(result.solution[1], 0.6, 1.0e-12, "scaled solve x1");
    }

    DenseMatrix huge{{1.0e308, 1.0e308}, {1.0e308, -1.0e308}};
    auto huge_result = solve_dense_partial_pivot(huge, Vector{1.0e308, 1.0});
    expect_true(huge_result.status == SolverStatus::breakdown || huge_result.status == SolverStatus::not_converged,
                "finite input numeric failure is structured");
}

void fixed_seed_random_cross_check() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::Vector;
    using MatCal::Linalg::factorize_dense_partial_pivot;
    using MatCal::Linalg::solve_dense_partial_pivot;

    std::mt19937 rng(7817);
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    for (int sample = 0; sample < 8; ++sample) {
        DenseMatrix a(4, 4);
        for (std::size_t r = 0; r < 4; ++r) {
            double row_sum = 0.0;
            for (std::size_t c = 0; c < 4; ++c) {
                if (r == c) {
                    continue;
                }
                const double value = dist(rng);
                a(r, c) = value;
                row_sum += std::abs(value);
            }
            a(r, r) = row_sum + 2.0 + 0.1 * static_cast<double>(sample);
        }
        Vector expected{1.0, -2.0, 0.5, 3.0};
        Vector rhs = a.multiply(expected);
        auto factorization = factorize_dense_partial_pivot(a);
        auto reused = factorization.factorization.solve(rhs);
        auto one_shot = solve_dense_partial_pivot(a, rhs);
        expect_true(factorization.success() && reused.success() && one_shot.success(), "random solve success");
        for (std::size_t i = 0; i < expected.size(); ++i) {
            expect_near(reused.solution[i], expected[i], 1.0e-11, "random reused solution");
            expect_near(reused.solution[i], one_shot.solution[i], 1.0e-12, "random reused vs one-shot");
        }
    }
}

} // namespace

int main() {
    factorization_contracts();
    determinant_contracts();
    failure_and_scale_contracts();
    fixed_seed_random_cross_check();
    return finish("MatCal::Linalg pivoted LU tests");
}
