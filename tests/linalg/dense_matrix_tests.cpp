#include <limits>

#include "MatCal/Linalg/DenseMatrix.hpp"
#include "linalg_test_support.hpp"

using namespace matcal_linalg_test;

int main() {
    using MatCal::Linalg::DenseMatrix;
    using MatCal::Linalg::Vector;

    DenseMatrix empty;
    expect_true(empty.rows() == 0 && empty.cols() == 0 && empty.empty(), "default matrix is 0x0");

    DenseMatrix zero_by_n(0, 4);
    expect_true(zero_by_n.rows() == 0 && zero_by_n.cols() == 4 && zero_by_n.size() == 0, "0xN matrix contract");

    DenseMatrix n_by_zero(3, 0);
    expect_true(n_by_zero.rows() == 3 && n_by_zero.cols() == 0 && n_by_zero.size() == 0, "Nx0 matrix contract");
    expect_true(n_by_zero.row(1).empty(), "Nx0 row view is empty");
    expect_throw([&] { (void)n_by_zero.row(3); }, "row view out of range");

    DenseMatrix a{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    expect_true(a.rows() == 2 && a.cols() == 3, "initializer dimensions");
    expect_near(a(1, 2), 6.0, 1e-12, "operator() access");
    a.at(1, 2) = 7.0;
    expect_near(a.at(1, 2), 7.0, 1e-12, "checked at access");
    expect_throw([&] { (void)a.at(2, 0); }, "DenseMatrix at row out of range");
    expect_throw([&] { DenseMatrix bad{{1.0}, {2.0, 3.0}}; }, "ragged initializer rejected");

    auto row = a.row(0);
    row[1] = 8.0;
    expect_near(a(0, 1), 8.0, 1e-12, "writable row view");
    std::span<const double> const_row = static_cast<const DenseMatrix&>(a).row(0);
    expect_true(const_row.size() == 3, "const row view");

    DenseMatrix identity = DenseMatrix::identity(3);
    expect_near(identity(0, 0), 1.0, 1e-12, "identity diagonal");
    expect_near(identity(0, 1), 0.0, 1e-12, "identity off diagonal");

    DenseMatrix transposed = a.transpose();
    expect_true(transposed.rows() == 3 && transposed.cols() == 2, "transpose dimensions");
    expect_near(transposed(2, 1), 7.0, 1e-12, "transpose data");

    Vector x{1.0, 2.0, 3.0};
    Vector ax = a.multiply(x);
    expect_near(ax[0], 26.0, 1e-12, "matvec row 0");
    expect_near(ax[1], 35.0, 1e-12, "matvec row 1");
    expect_throw([&] { (void)a.multiply(Vector{1.0, 2.0}); }, "matvec dimension mismatch");

    DenseMatrix b{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    DenseMatrix ab = a.multiply(b);
    expect_true(ab.rows() == 2 && ab.cols() == 2, "matmul dimensions");
    expect_near(ab(0, 0), 40.0, 1e-12, "matmul 00");
    expect_near(ab(1, 1), 70.0, 1e-12, "matmul 11");
    expect_throw([&] { (void)a.multiply(DenseMatrix{{1.0, 2.0}}); }, "matmul dimension mismatch");

    DenseMatrix filled(2, 2);
    filled.fill(3.0);
    expect_near(filled(1, 1), 3.0, 1e-12, "fill");

    expect_throw([] {
        DenseMatrix too_large(std::numeric_limits<std::size_t>::max() / 2 + 1, 3);
    }, "size overflow rejected before allocation");

    DenseMatrix non_finite{{1.0, std::numeric_limits<double>::infinity()}};
    expect_true(!non_finite.all_finite(), "finite check detects Inf");
    DenseMatrix nan_matrix{{std::numeric_limits<double>::quiet_NaN()}};
    expect_true(!nan_matrix.all_finite(), "finite check detects NaN");

    DenseMatrix copy = a;
    copy(0, 0) = 99.0;
    expect_near(a(0, 0), 1.0, 1e-12, "copy independence");
    DenseMatrix moved = std::move(copy);
    expect_near(moved(0, 0), 99.0, 1e-12, "move construction");

    return finish("MatCal::Linalg DenseMatrix tests");
}
