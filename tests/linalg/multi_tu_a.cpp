#include "MatCal/Linalg/DenseSolver.hpp"

MatCal::Linalg::Vector matcal_linalg_multi_tu_solve_a() {
    auto result = MatCal::Linalg::solve_dense_partial_pivot(
        MatCal::Linalg::DenseMatrix{{2.0, 0.0}, {0.0, 4.0}},
        MatCal::Linalg::Vector{6.0, 8.0});
    return result.solution;
}
