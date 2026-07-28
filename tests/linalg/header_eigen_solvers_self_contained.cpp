#include "MatCal/Linalg/EigenSolvers.hpp"

int header_eigen_solvers_self_contained() {
    auto result = MatCal::Linalg::dominant_eigenpair(MatCal::Linalg::DenseMatrix::identity(1));
    return result.success() ? 0 : 1;
}
