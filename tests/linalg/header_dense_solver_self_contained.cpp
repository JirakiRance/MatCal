#include "MatCal/Linalg/DenseSolver.hpp"

int matcal_linalg_header_dense_solver_probe() {
    auto result = MatCal::Linalg::solve_dense_partial_pivot(
        MatCal::Linalg::DenseMatrix::identity(1),
        MatCal::Linalg::Vector{2.0});
    return result.success() ? 0 : 1;
}
