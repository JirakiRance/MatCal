#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/SkylineLdlt.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

namespace {

MatCal::Linalg::SymmetricSkylineMatrix make_profile_spd(const std::vector<std::size_t>& first_columns,
                                                        double scale = 1.0) {
    using MatCal::Linalg::SymmetricSkylineMatrix;

    SymmetricSkylineMatrix matrix = SymmetricSkylineMatrix::from_profile(first_columns);
    std::vector<double> row_abs(first_columns.size(), 0.0);

    for (std::size_t row = 0; row < first_columns.size(); ++row) {
        for (std::size_t column = first_columns[row]; column < row; ++column) {
            double value = -scale * (0.01 + 0.002 * static_cast<double>((row + column) % 5));
            matrix.set(row, column, value);
            row_abs[row] += std::abs(value);
            row_abs[column] += std::abs(value);
        }
    }

    for (std::size_t row = 0; row < first_columns.size(); ++row) {
        matrix.set(row, row, scale * (1.0 + row_abs[row] / std::max(scale, 1.0)));
    }
    return matrix;
}

MatCal::Linalg::DenseMatrix to_dense(const MatCal::Linalg::SymmetricSkylineMatrix& matrix) {
    MatCal::Linalg::DenseMatrix dense(matrix.size(), matrix.size());
    for (std::size_t row = 0; row < matrix.size(); ++row) {
        for (std::size_t column = 0; column < matrix.size(); ++column) {
            dense(row, column) = matrix.get(row, column);
        }
    }
    return dense;
}

MatCal::Linalg::Vector deterministic_x(std::size_t n, double offset = 0.0) {
    MatCal::Linalg::Vector x(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = offset + 1.0 + static_cast<double>(i % 7) * 0.25;
    }
    return x;
}

std::vector<std::size_t> band_profile(std::size_t n, std::size_t half_band) {
    std::vector<std::size_t> first_columns(n);
    for (std::size_t row = 0; row < n; ++row) {
        first_columns[row] = row > half_band ? row - half_band : 0;
    }
    return first_columns;
}

std::vector<std::size_t> variable_profile(std::size_t n) {
    std::vector<std::size_t> first_columns(n);
    for (std::size_t row = 0; row < n; ++row) {
        std::size_t width = 1 + (row % 6);
        first_columns[row] = row + 1 > width ? row + 1 - width : 0;
    }
    return first_columns;
}

void expect_skyline_dense_agree(const MatCal::Linalg::SymmetricSkylineMatrix& skyline,
                                const MatCal::Linalg::Vector& rhs,
                                const std::string& label) {
    auto skyline_factor = MatCal::Linalg::factorize_skyline_ldlt(skyline);
    expect_true(skyline_factor.success(), label + " skyline factorization");
    auto skyline_result = skyline_factor.factorization.solve(rhs);
    auto dense_result = MatCal::Linalg::solve_dense_partial_pivot(to_dense(skyline), rhs);
    expect_true(skyline_result.success(), label + " skyline solve");
    expect_true(dense_result.success(), label + " dense solve");
    expect_true(skyline_result.metrics.factorization_operation_count ==
                    skyline_factor.metrics.factorization_operation_count,
                label + " solve reports existing factorization work");
    expect_true(skyline_result.metrics.solve_operation_count > 0, label + " solve operation count");
    expect_true(skyline_result.metrics.operation_count ==
                    skyline_result.metrics.factorization_operation_count +
                        skyline_result.metrics.solve_operation_count,
                label + " total operation count");
    expect_true(skyline_result.metrics.absolute_residual_norm <=
                    skyline_result.metrics.residual_acceptance_tolerance,
                label + " skyline residual accepted");
    expect_true(dense_result.metrics.absolute_residual_norm <= dense_result.metrics.residual_acceptance_tolerance,
                label + " dense residual accepted");
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        expect_near(skyline_result.solution[i], dense_result.solution[i], 1e-9, label + " solution");
    }
    double residual_tolerance = 1e-10 * std::max(skyline_result.metrics.matrix_scale, skyline_result.metrics.rhs_scale);
    residual_tolerance = std::max(residual_tolerance, 1e-12);
    expect_near(skyline_result.metrics.absolute_residual_norm,
                dense_result.metrics.absolute_residual_norm,
                residual_tolerance,
                label + " absolute residual");
    expect_near(skyline_result.metrics.relative_residual_norm,
                dense_result.metrics.relative_residual_norm,
                1e-10,
                label + " relative residual");
}

} // namespace

int main() {
    using MatCal::Linalg::SolverOptions;
    using MatCal::Linalg::SymmetricSkylineMatrix;

    auto narrow = make_profile_spd(band_profile(48, 2));
    auto varying = make_profile_spd(variable_profile(48));
    auto compact = make_profile_spd(band_profile(12, 1));

    expect_true(narrow.storage_size() < narrow.size() * narrow.size() / 8,
                "narrow skyline storage is far below dense");
    expect_true(narrow.storage_size() == 141, "narrow storage baseline");
    expect_true(varying.storage_size() > narrow.storage_size(), "variable profile stores more than narrow band");
    expect_true(varying.storage_size() == 168, "variable storage baseline");
    expect_true(varying.storage_size() < varying.size() * varying.size(), "variable profile stays below dense storage");

    auto rhs = narrow.multiply(deterministic_x(narrow.size()));
    auto factor = MatCal::Linalg::factorize_skyline_ldlt(narrow);
    expect_true(factor.success(), "narrow factorization succeeds");
    expect_true(factor.factorization.storage_size() == narrow.storage_size(), "factor storage follows skyline");
    expect_true(factor.metrics.factorization_operation_count == factor.metrics.operation_count,
                "factorization operation count recorded");
    expect_true(factor.metrics.factorization_operation_count == 139, "narrow factorization work baseline");
    expect_true(factor.factorization.factorization_operation_count() ==
                    factor.metrics.factorization_operation_count,
                "factorization owns operation count");

    auto solve_a = factor.factorization.solve(rhs);
    auto solve_b = factor.factorization.solve(narrow.multiply(deterministic_x(narrow.size(), -0.5)));
    expect_true(solve_a.success(), "first reused solve succeeds");
    expect_true(solve_b.success(), "second reused solve succeeds");
    expect_true(solve_a.metrics.factorization_operation_count == factor.metrics.factorization_operation_count,
                "first solve reuses factorization metrics");
    expect_true(solve_b.metrics.factorization_operation_count == factor.metrics.factorization_operation_count,
                "second solve reuses factorization metrics");
    expect_true(solve_a.metrics.solve_operation_count == solve_b.metrics.solve_operation_count,
                "same profile gives deterministic solve work");
    expect_true(solve_a.metrics.solve_operation_count == 375, "narrow solve work uses skyline profile");

    expect_skyline_dense_agree(compact, compact.multiply(deterministic_x(compact.size())), "compact exact");
    expect_skyline_dense_agree(varying, varying.multiply(deterministic_x(varying.size())), "varying profile");

    for (double scale : {1e-20, 1.0, 1e20}) {
        auto scaled = make_profile_spd(band_profile(8, 2), scale);
        expect_skyline_dense_agree(scaled, scaled.multiply(deterministic_x(scaled.size())), "scaled " + std::to_string(scale));
    }

    SolverOptions caller_options;
    std::size_t n = compact.size();
    double c = 256.0 * std::numeric_limits<double>::epsilon() * static_cast<double>(n);
    caller_options.absolute_tolerance = c;
    caller_options.relative_tolerance = c;
    caller_options.pivot_factor = 1.0;
    auto caller_factor = MatCal::Linalg::factorize_skyline_ldlt(compact, caller_options);
    expect_true(caller_factor.success(), "caller tolerance expression factorizes");
    double expected_tolerance = c * std::max(caller_factor.metrics.matrix_scale, 1.0);
    expect_near(caller_factor.metrics.pivot_tolerance_used,
                expected_tolerance,
                expected_tolerance * 1e-12,
                "caller can express c times max(matrix_scale, 1)");

    SymmetricSkylineMatrix one(1);
    one.set(0, 0, 3.0);
    auto one_factor = MatCal::Linalg::factorize_skyline_ldlt(one);
    expect_true(one_factor.success(), "1x1 performance factorization");
    expect_true(one_factor.metrics.factorization_operation_count == 0, "1x1 factorization has no offdiagonal work");
    auto one_solve = one_factor.factorization.solve(MatCal::Linalg::Vector{6.0});
    expect_true(one_solve.success(), "1x1 performance solve");
    expect_true(one_solve.metrics.solve_operation_count > 0, "1x1 solve work counted");

    return finish("MatCal::Linalg Skyline performance contract tests");
}
