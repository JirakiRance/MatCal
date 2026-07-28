#include <cmath>
#include <limits>

#include "test_support.hpp"
#include "Basics.hpp"
#include "Matrix.hpp"

namespace {

using matcal_test::expect_near;
using matcal_test::expect_throw;
using matcal_test::expect_true;
using matcal_test::finish;

void conversion_adapters_deep_copy() {
    MatCal::Utils::Matrix legacy{{1.0, 2.0}, {3.0, 4.0}};
    auto dense = MatCal::Legacy::to_linalg_dense(legacy);
    legacy.set(0, 0, 9.0);
    expect_near(dense(0, 0), 1.0, 0.0, "legacy to linalg is deep copy");

    dense(0, 1) = 8.0;
    auto legacy_copy = MatCal::Legacy::to_legacy_matrix(dense);
    dense(0, 1) = 7.0;
    expect_near(legacy_copy.get(0, 1), 8.0, 0.0, "linalg to legacy is deep copy");

    MatCal::Utils::Matrix column{{5.0}, {6.0}};
    auto vector = MatCal::Legacy::to_linalg_column_vector(column);
    auto back = MatCal::Legacy::to_legacy_column_matrix(vector);
    expect_near(back.get(1, 0), 6.0, 0.0, "column vector roundtrip");
    expect_throw([&] { (void)MatCal::Legacy::to_linalg_column_vector(legacy); },
                 "adapter rejects non-column vector");

    MatCal::Utils::Matrix empty(0, 0);
    auto empty_dense = MatCal::Legacy::to_linalg_dense(empty);
    auto empty_back = MatCal::Legacy::to_legacy_matrix(empty_dense);
    expect_true(empty_dense.rows() == 0 && empty_dense.cols() == 0, "empty matrix converts to linalg");
    expect_true(empty_back.getRows() == 0 && empty_back.getCols() == 0, "empty matrix roundtrip");

    MatCal::Utils::Matrix nonfinite{{1.0, std::numeric_limits<double>::quiet_NaN()}};
    expect_throw([&] { (void)MatCal::Legacy::to_linalg_dense(nonfinite); },
                 "adapter rejects non-finite matrix values");
}

void direct_solver_facade_matches_linalg() {
    MatCal::Utils::Matrix a{{0.0, 1.0}, {1.0, 1.0}};
    MatCal::Utils::Matrix b{{2.0}, {5.0}};
    auto solved = MatCal::Algorithm::Matrix::solve_columnElimination(a, b);
    expect_near(solved->get(0, 0), 3.0, 1.0e-12, "legacy direct solver x0");
    expect_near(solved->get(1, 0), 2.0, 1.0e-12, "legacy direct solver x1");
    expect_near(a.get(0, 0), 0.0, 0.0, "legacy direct solver does not modify A");

    MatCal::Utils::Matrix multi_rhs{{2.0, 4.0}, {5.0, 11.0}};
    auto multi = MatCal::Algorithm::Matrix::solve_columnElimination(a, multi_rhs);
    expect_near(multi->get(0, 0), 3.0, 1.0e-12, "legacy direct multi RHS col0 x0");
    expect_near(multi->get(1, 0), 2.0, 1.0e-12, "legacy direct multi RHS col0 x1");
    expect_near(multi->get(0, 1), 7.0, 1.0e-12, "legacy direct multi RHS col1 x0");
    expect_near(multi->get(1, 1), 4.0, 1.0e-12, "legacy direct multi RHS col1 x1");

    MatCal::Utils::Matrix empty_a(0, 0);
    MatCal::Utils::Matrix empty_b(0, 0);
    auto empty_solution = MatCal::Algorithm::Matrix::solve_columnElimination(empty_a, empty_b);
    expect_true(empty_solution->getRows() == 0 && empty_solution->getCols() == 0,
                "legacy direct empty system returns empty matrix");

    MatCal::Utils::Matrix rhs_with_nan{{1.0, std::numeric_limits<double>::quiet_NaN()}, {2.0, 3.0}};
    expect_throw([&] { (void)MatCal::Algorithm::Matrix::solve_columnElimination(a, rhs_with_nan); },
                 "legacy direct rejects non-finite RHS before returning a partial result");
    expect_near(a.get(0, 0), 0.0, 0.0, "legacy direct failed multi RHS keeps A unchanged");
    expect_near(rhs_with_nan.get(1, 1), 3.0, 0.0, "legacy direct failed multi RHS keeps b unchanged");
}

void determinant_and_lu_compatibility_edges() {
    MatCal::Utils::Matrix swap_matrix{{0.0, 1.0}, {1.0, 0.0}};
    double det = MatCal::Algorithm::Matrix::determinant(swap_matrix);
    expect_near(det, -1.0, 0.0, "determinant keeps row-swap parity");
    expect_near(swap_matrix.get(0, 0), 0.0, 0.0, "determinant does not modify input");

    MatCal::Utils::Matrix pivoting_required{{0.0, 1.0}, {1.0, 1.0}};
    expect_throw([&] { (void)MatCal::Algorithm::Matrix::LU_Decompose(pivoting_required); },
                 "no-pivot LU reports pivoting-required failure");
}

void stationary_legacy_delegates_linalg() {
    MatCal::Utils::Matrix a{{4.0, 1.0}, {2.0, 3.0}};
    MatCal::Utils::Matrix b{{1.0}, {2.0}};
    auto dense = MatCal::Legacy::to_linalg_dense(a);
    auto rhs = MatCal::Legacy::to_linalg_column_vector(b);
    MatCal::Linalg::SolverOptions options;
    options.absolute_tolerance = 0.0;
    options.relative_tolerance = 1.0e-8;
    options.max_iterations = 200;

    auto core = MatCal::Linalg::solve_gauss_seidel(dense, rhs, options);
    auto legacy = MatCal::Algorithm::Matrix::Gauss_Seidel(a, b, 1.0e-8, 200);
    expect_true(core.success(), "core GS success");
    expect_true(legacy.converged, "legacy GS success");
    expect_near(legacy.root.get(0, 0), core.solution[0], 1.0e-8, "legacy/core GS x0");
    expect_near(legacy.root.get(1, 0), core.solution[1], 1.0e-8, "legacy/core GS x1");

    auto sor = MatCal::Algorithm::Matrix::SOR(a, b, 1.2, 1.0e-8, 200);
    expect_true(sor.converged, "legacy SOR double delegates");
    auto sor_int = MatCal::Algorithm::Matrix::SOR(a, b, 1, 1.0e-8, 200);
    expect_true(sor_int.converged, "legacy SOR int remains available");
}

void eigen_legacy_delegates_linalg() {
    MatCal::Utils::Matrix a{{2.0, 1.0}, {1.0, 2.0}};
    auto legacy = MatCal::Algorithm::Matrix::PowerMethod(a, 1.0e-8, 200);
    auto core = MatCal::Linalg::dominant_eigenpair(MatCal::Legacy::to_linalg_dense(a));
    expect_true(core.success(), "core power success");
    expect_near(legacy.eigenvalue, core.eigenvalue, 1.0e-8, "legacy/core power eigenvalue");

    MatCal::Utils::Matrix d{{3.0, 0.0}, {0.0, 2.0}};
    auto inverse = MatCal::Algorithm::Matrix::PowerMethod_reverse(d, 1.9, 1.0e-8, 200);
    expect_near(inverse.eigenvalue, 2.0, 1.0e-8, "legacy inverse power eigenvalue");
}

void derivative_helpers_delegate_calculus() {
    using Derivative = MatCal::Algorithm::Basics::Derivative;
    Derivative::Func_F f = [](std::vector<double> values) {
        return values[0] * values[0] + 3.0 * values[1];
    };
    std::vector<double> point{2.0, 4.0};
    expect_near(Derivative::pF_px(f, point, 0, 1.0e-6), 4.0, 1.0e-5, "legacy partial x0");
    expect_near(Derivative::pF_px(f, point, 1, 1.0e-6), 3.0, 1.0e-8, "legacy partial x1");
    expect_near(Derivative::dF_dx(f, point, 1.0e-6), 7.0, 1.0e-5, "legacy gradient sum");
    expect_near(point[0], 2.0, 0.0, "partial derivative does not mutate input");
    expect_throw([&] { (void)Derivative::pF_px(f, point, -1); }, "partial derivative rejects negative index");
}

} // namespace

int main() {
    conversion_adapters_deep_copy();
    direct_solver_facade_matches_linalg();
    determinant_and_lu_compatibility_edges();
    stationary_legacy_delegates_linalg();
    eigen_legacy_delegates_linalg();
    derivative_helpers_delegate_calculus();
    return finish("legacy Matrix/Linalg adapter tests");
}
