#include <limits>
#include <string>

#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/SkylineLdlt.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

MatCal::Linalg::SymmetricSkylineMatrix make_spd(double scale = 1.0) {
    auto matrix = MatCal::Linalg::SymmetricSkylineMatrix::from_profile({0, 0, 1, 2});
    matrix.set(0, 0, 4.0 * scale);
    matrix.set(1, 0, 1.0 * scale);
    matrix.set(1, 1, 3.0 * scale);
    matrix.set(2, 1, 0.5 * scale);
    matrix.set(2, 2, 2.0 * scale);
    matrix.set(3, 2, 1.0 * scale);
    matrix.set(3, 3, 2.0 * scale);
    return matrix;
}

MatCal::Linalg::DenseMatrix make_spd_dense(double scale = 1.0) {
    return MatCal::Linalg::DenseMatrix{{4.0 * scale, 1.0 * scale, 0.0, 0.0},
                                       {1.0 * scale, 3.0 * scale, 0.5 * scale, 0.0},
                                       {0.0, 0.5 * scale, 2.0 * scale, 1.0 * scale},
                                       {0.0, 0.0, 1.0 * scale, 2.0 * scale}};
}

void expect_known_solution(const MatCal::Linalg::SolverResult& result, const std::string& label) {
    expect_true(result.success(), label + " success");
    expect_near(result.solution[0], 1.0, 1e-10, label + " x0");
    expect_near(result.solution[1], 2.0, 1e-10, label + " x1");
    expect_near(result.solution[2], 3.0, 1e-10, label + " x2");
    expect_near(result.solution[3], 4.0, 1e-10, label + " x3");
    expect_true(result.metrics.absolute_residual_norm <= result.metrics.residual_acceptance_tolerance,
                label + " residual accepted");
}

} // namespace

int main() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::SolverOptions;
    using MatCal::Linalg::SolverStatus;
    using MatCal::Linalg::SymmetricSkylineMatrix;
    using MatCal::Linalg::Vector;
    using MatCal::Linalg::factorize_skyline_ldlt;
    using MatCal::Linalg::solve_dense_partial_pivot;
    using MatCal::Linalg::solve_skyline_ldlt;

    SymmetricSkylineMatrix empty;
    auto empty_factor = factorize_skyline_ldlt(empty);
    expect_true(empty_factor.success(), "0x0 factorization succeeds");
    auto empty_solve = empty_factor.factorization.solve(Vector());
    expect_true(empty_solve.success() && empty_solve.solution.empty(), "0x0 solve succeeds");

    SymmetricSkylineMatrix one(1);
    one.set(0, 0, 4.0);
    auto one_result = solve_skyline_ldlt(one, Vector{8.0});
    expect_true(one_result.success(), "1x1 solve succeeds");
    expect_near(one_result.solution[0], 2.0, 1e-12, "1x1 solution");

    auto matrix = make_spd();
    Vector x_expected{1.0, 2.0, 3.0, 4.0};
    Vector rhs = matrix.multiply(x_expected);
    auto original_10 = matrix.get(1, 0);
    auto factorized = factorize_skyline_ldlt(matrix);
    expect_true(factorized.success(), "small SPD factorization succeeds");
    expect_true(factorized.factorization.storage_size() == matrix.storage_size(), "factor owns skyline-sized storage");
    auto result = factorized.factorization.solve(rhs);
    expect_known_solution(result, "small SPD solve");
    expect_near(matrix.get(1, 0), original_10, 1e-12, "factorization does not modify input");

    auto dense_result = solve_dense_partial_pivot(make_spd_dense(), rhs);
    expect_true(dense_result.success(), "dense cross-check succeeds");
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        expect_near(result.solution[i], dense_result.solution[i], 1e-10, "skyline dense cross-check");
    }

    Vector rhs2 = matrix.multiply(Vector{2.0, -1.0, 0.5, 3.0});
    auto result2 = factorized.factorization.solve(rhs2);
    expect_true(result2.success(), "factorization reused for second RHS");
    expect_near(result2.solution[0], 2.0, 1e-10, "second RHS x0");
    expect_near(result2.solution[1], -1.0, 1e-10, "second RHS x1");
    expect_near(result2.solution[2], 0.5, 1e-10, "second RHS x2");
    expect_near(result2.solution[3], 3.0, 1e-10, "second RHS x3");

    for (double scale : {1e-20, 1.0, 1e20}) {
        auto scaled_matrix = make_spd(scale);
        auto scaled_rhs = scaled_matrix.multiply(x_expected);
        auto scaled_result = solve_skyline_ldlt(scaled_matrix, scaled_rhs);
        expect_known_solution(scaled_result, "scaled skyline " + std::to_string(scale));
    }

    SymmetricSkylineMatrix singular = SymmetricSkylineMatrix::from_profile({0, 0});
    singular.set(0, 0, 1.0);
    singular.set(1, 0, 1.0);
    singular.set(1, 1, 1.0);
    auto singular_result = factorize_skyline_ldlt(singular);
    expect_true(singular_result.status == SolverStatus::not_positive_definite, "singular rejected as non-SPD");
    expect_true(!singular_result.diagnostics.empty() && singular_result.diagnostics[0].code == "non_positive_pivot",
                "singular diagnostic code");

    SymmetricSkylineMatrix indefinite = SymmetricSkylineMatrix::from_profile({0, 0});
    indefinite.set(0, 0, 1.0);
    indefinite.set(1, 0, 2.0);
    indefinite.set(1, 1, 1.0);
    auto indefinite_result = factorize_skyline_ldlt(indefinite);
    expect_true(indefinite_result.status == SolverStatus::not_positive_definite, "indefinite rejected");

    SymmetricSkylineMatrix non_finite(1);
    non_finite.set(0, 0, std::numeric_limits<double>::quiet_NaN());
    auto non_finite_matrix = factorize_skyline_ldlt(non_finite);
    expect_true(non_finite_matrix.status == SolverStatus::non_finite_input, "non-finite matrix rejected");
    auto non_finite_rhs = factorized.factorization.solve(Vector{1.0, 2.0, std::numeric_limits<double>::infinity(), 4.0});
    expect_true(non_finite_rhs.status == SolverStatus::non_finite_input, "non-finite RHS rejected");

    SymmetricSkylineMatrix overflow = SymmetricSkylineMatrix::from_profile({0, 0});
    overflow.set(0, 0, 1e-100);
    overflow.set(1, 0, 1e200);
    overflow.set(1, 1, 1.0);
    SolverOptions exact_pivot;
    exact_pivot.absolute_tolerance = 0.0;
    exact_pivot.relative_tolerance = 0.0;
    auto overflow_factor = factorize_skyline_ldlt(overflow, exact_pivot);
    expect_true(overflow_factor.status == SolverStatus::breakdown, "finite factorization overflow returns breakdown");
    expect_true(!overflow_factor.diagnostics.empty() && overflow_factor.diagnostics[0].code == "non_finite_intermediate",
                "overflow diagnostic code");

    SymmetricSkylineMatrix gap = SymmetricSkylineMatrix::from_profile({0, 0, 0});
    gap.set(0, 0, 2.0);
    gap.set(1, 1, 2.0);
    gap.set(2, 0, 0.25);
    gap.set(2, 2, 2.0);
    expect_near(gap.get(2, 1), 0.0, 1e-12, "profile gap defaults to zero");
    auto gap_result = solve_skyline_ldlt(gap, gap.multiply(Vector{1.0, 2.0, 3.0}));
    expect_true(gap_result.success(), "profile gap matrix solves");

    SymmetricSkylineMatrix diagonal(4);
    expect_true(diagonal.storage_size() < 16, "skyline storage smaller than dense");

    SolverOptions strict;
    strict.absolute_tolerance = 0.0;
    strict.relative_tolerance = 0.0;
    auto strict_result = solve_skyline_ldlt(matrix, rhs, strict);
    expect_true(strict_result.success(), "strict exact tolerance solve succeeds for small SPD");
    expect_true(strict_result.metrics.matrix_scale > 0.0, "matrix scale recorded");
    expect_true(strict_result.metrics.pivot_tolerance_used == 0.0, "strict pivot tolerance recorded");
    expect_true(strict_result.metrics.minimum_abs_pivot > 0.0, "minimum pivot recorded");

    auto first = factorize_skyline_ldlt(singular);
    auto second = factorize_skyline_ldlt(singular);
    expect_true(first.status == second.status, "deterministic factorization status");
    expect_true(first.diagnostics[0].code == second.diagnostics[0].code, "deterministic diagnostic code");
    expect_true(first.diagnostics[0].phase == second.diagnostics[0].phase, "deterministic diagnostic phase");

    return finish("MatCal::Linalg Skyline LDLT tests");
}
