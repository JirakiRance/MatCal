#include "MatCal/Linalg/SkylineLdlt.hpp"

double matcal_linalg_multi_tu_skyline_value() {
    MatCal::Linalg::SymmetricSkylineMatrix matrix(1);
    matrix.set(0, 0, 5.0);
    auto result = MatCal::Linalg::solve_skyline_ldlt(matrix, MatCal::Linalg::Vector{10.0});
    return result.solution[0];
}
