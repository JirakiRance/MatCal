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
}

void direct_solver_facade_matches_linalg() {
    MatCal::Utils::Matrix a{{0.0, 1.0}, {1.0, 1.0}};
    MatCal::Utils::Matrix b{{2.0}, {5.0}};
    auto solved = MatCal::Algorithm::Matrix::solve_columnElimination(a, b);
    expect_near(solved->get(0, 0), 3.0, 1.0e-12, "legacy direct solver x0");
    expect_near(solved->get(1, 0), 2.0, 1.0e-12, "legacy direct solver x1");
    expect_near(a.get(0, 0), 0.0, 0.0, "legacy direct solver does not modify A");
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
    stationary_legacy_delegates_linalg();
    eigen_legacy_delegates_linalg();
    derivative_helpers_delegate_calculus();
    return finish("legacy Matrix/Linalg adapter tests");
}
