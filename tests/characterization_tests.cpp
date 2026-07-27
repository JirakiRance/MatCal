#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Basics.hpp"
#include "Insert.hpp"
#include "Iteration.hpp"
#include "Matrix.hpp"
#include "QinJiuShao.hpp"

namespace {

int failures = 0;

void fail(const std::string& message) {
    ++failures;
    std::cerr << "FAIL: " << message << '\n';
}

void expect_true(bool condition, const std::string& message) {
    if (!condition) {
        fail(message);
    }
}

void expect_near(double actual, double expected, double tolerance, const std::string& message) {
    if (std::isnan(expected)) {
        if (!std::isnan(actual)) {
            fail(message + " expected NaN");
        }
        return;
    }
    if (std::abs(actual - expected) > tolerance) {
        fail(message + " expected " + std::to_string(expected) + " got " + std::to_string(actual));
    }
}

template <typename F>
void expect_throw(F&& func, const std::string& message) {
    try {
        func();
    } catch (const std::exception&) {
        return;
    }
    fail(message + " expected exception");
}

MatCal::Utils::Matrix as_matrix(std::unique_ptr<MatCal::Utils::AbstractMatrix> value) {
    auto* matrix = dynamic_cast<MatCal::Utils::Matrix*>(value.get());
    if (!matrix) {
        throw std::runtime_error("Expected legacy operation to return MatCal::Utils::Matrix");
    }
    return *matrix;
}

void matrix_basics_are_characterized() {
    using MatCal::Utils::Matrix;

    Matrix raw({{1.0, 2.0}, {3.0, 4.0}});
    expect_true(raw.getRows() == 2 && raw.getCols() == 2, "Matrix initializer-list dimensions");
    expect_near(raw.get(1, 0), 3.0, 1e-12, "Matrix get");
    raw.set(0, 1, 5.0);
    expect_near(raw[0][1], 5.0, 1e-12, "Matrix set and operator[]");

    const double r0[] = {1.0, 2.0};
    const double r1[] = {3.0, 4.0};
    const double* rows[] = {r0, r1};
    const double** ptr = rows;
    Matrix from_raw(ptr, 2, 2);
    expect_near(from_raw.get(1, 1), 4.0, 1e-12, "Matrix raw double** constructor copies values");

    Matrix a({{1.0, 2.0}, {3.0, 4.0}});
    Matrix b({{5.0, 6.0}, {7.0, 8.0}});
    Matrix sum = as_matrix(a.add(b));
    expect_near(sum.get(0, 0), 6.0, 1e-12, "Matrix addition");
    expect_near(sum.get(1, 1), 12.0, 1e-12, "Matrix addition lower-right");

    Matrix product = as_matrix(a.multiply(b));
    expect_near(product.get(0, 0), 19.0, 1e-12, "Matrix multiplication");
    expect_near(product.get(1, 1), 50.0, 1e-12, "Matrix multiplication lower-right");

    Matrix transposed = as_matrix(a.transpose());
    expect_near(transposed.get(0, 1), 3.0, 1e-12, "Matrix transpose");

    expect_throw([] { MatCal::Utils::Matrix invalid(-1, 2); }, "Negative Matrix dimensions");
    expect_throw([&] { (void)a.get(2, 0); }, "Matrix get out of bounds");
    expect_throw([&] { (void)a[2]; }, "Matrix operator[] out of bounds");
    expect_throw([&] {
        Matrix c(3, 1);
        (void)a.add(c);
    }, "Matrix add incompatible dimensions");
    expect_throw([] {
        const double** null_input = nullptr;
        MatCal::Utils::Matrix invalid(null_input, 2, 2);
    }, "Matrix null raw input");
}

void direct_solvers_are_characterized() {
    using namespace MatCal::Algorithm::Matrix;
    using MatCal::Utils::Matrix;
    using MatCal::Utils::TridiagonalMatrix;

    Matrix a({{2.0, 1.0}, {5.0, 7.0}});
    Matrix b({{11.0}, {13.0}});
    Matrix x = as_matrix(solve_columnElimination(a, b));
    expect_near(x.get(0, 0), 64.0 / 9.0, 1e-10, "Gaussian elimination solution x0");
    expect_near(x.get(1, 0), -29.0 / 9.0, 1e-10, "Gaussian elimination solution x1");
    expect_near(a.get(0, 0), 2.0, 1e-12, "Gaussian elimination does not modify A");
    expect_near(b.get(0, 0), 11.0, 1e-12, "Gaussian elimination does not modify b");

    Matrix singular({{1.0, 2.0}, {2.0, 4.0}});
    Matrix singular_rhs({{3.0}, {6.0}});
    expect_throw([&] { (void)solve_columnElimination(singular, singular_rhs); }, "Gaussian elimination singular matrix");

    Matrix lu_a({{4.0, 3.0}, {6.0, 3.0}});
    Matrix lu_b({{10.0}, {12.0}});
    auto lu = LU_Decompose(lu_a);
    Matrix lu_x = as_matrix(lu.solve(lu_b));
    expect_near(lu_x.get(0, 0), 1.0, 1e-10, "LU solve x0");
    expect_near(lu_x.get(1, 0), 2.0, 1e-10, "LU solve x1");

    Matrix pivot_needed({{0.0, 1.0}, {1.0, 1.0}});
    expect_throw([&] { (void)LU_Decompose(pivot_needed); }, "LU without pivoting fails on invertible zero-pivot matrix");

    TridiagonalMatrix tri({1.0, 1.0}, {4.0, 4.0, 4.0}, {1.0, 1.0});
    Matrix tri_rhs({{6.0}, {12.0}, {14.0}});
    Matrix tri_x = as_matrix(tri.solve(tri_rhs));
    expect_near(tri_x.get(0, 0), 1.0, 1e-10, "Thomas solve x0");
    expect_near(tri_x.get(1, 0), 2.0, 1e-10, "Thomas solve x1");
    expect_near(tri_x.get(2, 0), 3.0, 1e-10, "Thomas solve x2");
}

void iterative_solvers_are_characterized() {
    using namespace MatCal::Algorithm::Matrix;
    using MatCal::Utils::Matrix;

    Matrix a({{4.0, 1.0}, {2.0, 3.0}});
    Matrix b({{1.0}, {2.0}});

    auto jacobi = Jacobi(a, b, 1e-8, 200);
    expect_true(jacobi.converged, "Jacobi converges for diagonally dominant legacy case");
    expect_near(jacobi.root.get(0, 0), 0.1, 1e-5, "Jacobi x0");
    expect_near(jacobi.root.get(1, 0), 0.6, 1e-5, "Jacobi x1");

    auto gs = Gauss_Seidel(a, b, 1e-8, 200);
    expect_true(gs.converged, "Gauss-Seidel converges for diagonally dominant legacy case");
    expect_near(gs.root.get(0, 0), 0.1, 1e-5, "Gauss-Seidel x0");
    expect_near(gs.root.get(1, 0), 0.6, 1e-5, "Gauss-Seidel x1");

    auto sor = SOR(a, b, 1, 1e-8, 200);
    expect_true(sor.converged, "SOR with legacy int omega=1 converges");
    expect_near(sor.root.get(0, 0), 0.1, 1e-5, "SOR x0 with omega=1");
    expect_near(sor.root.get(1, 0), 0.6, 1e-5, "SOR x1 with omega=1");

    expect_throw([&] { (void)SOR(a, b, 2); }, "SOR rejects omega outside (0, 2)");
}

void polynomial_and_interpolation_are_characterized() {
    using MatCal::Algorithm::Insert::LagrangeInsert;
    using MatCal::Algorithm::Insert::LinearInsert;
    using MatCal::Utils::QinJiuShao;

    QinJiuShao poly({{2, 3.0}, {1, 2.0}, {0, 1.0}});
    expect_near(poly.calculate(2.0), 17.0, 1e-12, "QinJiuShao polynomial evaluation");
    expect_near(poly.derivative().calculate(2.0), 14.0, 1e-12, "QinJiuShao derivative");
    auto func = poly.toFunction();
    expect_near(func(2.0), 17.0, 1e-12, "QinJiuShao toFunction while owner is alive");

    LagrangeInsert lagrange({{0.0, 1.0}, {1.0, 3.0}, {2.0, 7.0}});
    expect_near(lagrange.calculate(3.0), 13.0, 1e-10, "Lagrange interpolation");

    LinearInsert linear({{0.0, 0.0}, {2.0, 4.0}});
    expect_near(linear.calculate(1.0), 2.0, 1e-12, "Linear interpolation");
    expect_near(linear.calculate(3.0), 6.0, 1e-12, "Linear extrapolation legacy behavior");
}

void nonlinear_and_basic_numerics_are_characterized() {
    using MatCal::Algorithm::Basics::NumericalIntegration;
    using MatCal::Algorithm::Iteration::Bisection;
    using MatCal::Algorithm::Iteration::Newton;

    auto root = Bisection::solveDetailed([](double x) { return x * x - 2.0; }, 1.0, 2.0, 1e-8, 200);
    expect_true(root.converged, "Bisection converges");
    expect_near(root.root, std::sqrt(2.0), 1e-7, "Bisection root");

    auto newton = Newton::solve([](double x) { return x * x - 2.0; },
                                [](double x) { return 2.0 * x; },
                                1.0, 1e-8, 200);
    expect_true(newton.converged, "Newton converges");
    expect_near(newton.root, std::sqrt(2.0), 1e-7, "Newton root");

    double integral = NumericalIntegration::NewtonCotes([](double x) { return x * x; }, 0.0, 1.0, 2);
    expect_near(integral, 1.0 / 3.0, 1e-12, "Newton-Cotes Simpson integration");
}

void nan_inf_behavior_is_characterized() {
    using MatCal::Utils::Matrix;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    Matrix nan_matrix({{nan}});
    Matrix nan_scaled = as_matrix(nan_matrix.scalarMultiply(2.0));
    expect_true(std::isnan(nan_scaled.get(0, 0)), "NaN propagates through scalarMultiply");

    Matrix inf_matrix({{inf}});
    Matrix inf_added = as_matrix(inf_matrix.add(Matrix({{1.0}})));
    expect_true(std::isinf(inf_added.get(0, 0)), "Inf propagates through add");
}

} // namespace

int main() {
    matrix_basics_are_characterized();
    direct_solvers_are_characterized();
    iterative_solvers_are_characterized();
    polynomial_and_interpolation_are_characterized();
    nonlinear_and_basic_numerics_are_characterized();
    nan_inf_behavior_is_characterized();

    if (failures != 0) {
        std::cerr << failures << " characterization test(s) failed\n";
        return 1;
    }
    std::cout << "MatCal legacy characterization tests passed\n";
    return 0;
}
