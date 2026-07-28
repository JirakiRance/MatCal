#include "MatCal/Linalg/IterativeSolvers.hpp"

int header_iterative_solvers_self_contained() {
    auto result = MatCal::Linalg::solve_jacobi(MatCal::Linalg::DenseMatrix::identity(1),
                                              MatCal::Linalg::Vector{1.0});
    return result.success() ? 0 : 1;
}
