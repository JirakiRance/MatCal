#include "MatCal/Linalg/DenseMatrix.hpp"

int matcal_linalg_header_dense_matrix_probe() {
    MatCal::Linalg::DenseMatrix m(1, 1);
    m(0, 0) = 1.0;
    return m.rows() == 1 ? 0 : 1;
}
