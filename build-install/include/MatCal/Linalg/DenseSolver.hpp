#ifndef MATCAL_LINALG_DENSE_SOLVER_HPP
#define MATCAL_LINALG_DENSE_SOLVER_HPP

#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/SolverTypes.hpp"
#include "MatCal/Linalg/Vector.hpp"

namespace MatCal::Linalg {

double residual_norm_inf(const DenseMatrix& a, const Vector& x, const Vector& b);

SolverResult solve_dense_partial_pivot(const DenseMatrix& a,
                                       const Vector& b,
                                       const SolverOptions& options = {});

} // namespace MatCal::Linalg

#endif
