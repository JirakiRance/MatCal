#include "MatCal/Linalg/DenseMatrix.hpp"

double matcal_linalg_multi_tu_matrix_value_b() {
    MatCal::Linalg::DenseMatrix matrix{{1.0, 2.0}, {3.0, 4.0}};
    return matrix.transpose()(0, 1);
}
