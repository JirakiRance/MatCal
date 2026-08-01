#include "MatCal/Linalg/DenseSolver.hpp"

int matcal_linalg_header_dense_solver_probe() {
    auto result = MatCal::Linalg::solve_dense_partial_pivot(
        MatCal::Linalg::DenseMatrix::identity(1),
        MatCal::Linalg::Vector{2.0});
    auto factorization = MatCal::Linalg::factorize_dense_partial_pivot(
        MatCal::Linalg::DenseMatrix{{0.0, 1.0}, {1.0, 1.0}});
    auto multi = factorization.factorization.solve(
        MatCal::Linalg::DenseMatrix{{1.0, 2.0}, {2.0, 3.0}});
    return result.success() && factorization.success() && multi.success() ? 0 : 1;
}
