#include "MatCal/Linalg/SkylineLdlt.hpp"

int matcal_linalg_header_skyline_ldlt_probe() {
    MatCal::Linalg::SymmetricSkylineMatrix matrix(1);
    matrix.set(0, 0, 2.0);
    auto result = MatCal::Linalg::solve_skyline_ldlt(matrix, MatCal::Linalg::Vector{4.0});
    return result.success() ? 0 : 1;
}
