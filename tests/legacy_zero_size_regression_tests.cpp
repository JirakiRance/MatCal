#include "test_support.hpp"

using namespace matcal_test;

int main() {
    using MatCal::Utils::Matrix;
    using MatCal::Utils::TridiagonalMatrix;

    TridiagonalMatrix zero(0);
    expect_true(zero.getRows() == 0 && zero.getCols() == 0, "zero tridiagonal dimensions");
    expect_true(zero.lower().empty(), "zero tridiagonal lower empty");
    expect_true(zero.diag().empty(), "zero tridiagonal diag empty");
    expect_true(zero.upper().empty(), "zero tridiagonal upper empty");

    Matrix dense = as_matrix(zero.toNormalMatrix());
    expect_true(dense.getRows() == 0 && dense.getCols() == 0, "zero tridiagonal dense conversion");

    Matrix rhs(0, 1);
    Matrix solved = as_matrix(zero.solve(rhs));
    expect_true(solved.getRows() == 0 && solved.getCols() == 1, "zero tridiagonal solve returns empty solution");

    Matrix other(0, 2);
    Matrix product = as_matrix(zero.multiply(other));
    expect_true(product.getRows() == 0 && product.getCols() == 2, "zero tridiagonal multiply is safe");

    expect_throw([&] { (void)zero.get(0, 0); }, "zero tridiagonal get out of range");
    return finish("legacy zero-size regression tests");
}
