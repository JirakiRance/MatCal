#include "MatCal/Linalg/EigenSolvers.hpp"
#include "MatCal/Linalg/IterativeSolvers.hpp"

int linalg_multi_tu_iterative_eigen() {
    auto linear = MatCal::Linalg::solve_gauss_seidel(MatCal::Linalg::DenseMatrix::identity(1),
                                                    MatCal::Linalg::Vector{2.0});
    auto eigen = MatCal::Linalg::dominant_eigenpair(MatCal::Linalg::DenseMatrix::identity(1));
    return linear.success() && eigen.success() ? 0 : 1;
}
