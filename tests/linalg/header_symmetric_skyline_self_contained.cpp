#include "MatCal/Linalg/SymmetricSkylineMatrix.hpp"

int matcal_linalg_header_symmetric_skyline_probe() {
    MatCal::Linalg::SymmetricSkylineMatrix matrix(1);
    matrix.set(0, 0, 1.0);
    return matrix.get(0, 0) == 1.0 ? 0 : 1;
}
