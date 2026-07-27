#include "MatCal/Linalg/Vector.hpp"

MatCal::Linalg::Vector matcal_linalg_multi_tu_solve_a();
double matcal_linalg_multi_tu_matrix_value_b();

int main() {
    auto solution = matcal_linalg_multi_tu_solve_a();
    if (solution.size() != 2 || solution[0] != 3.0 || solution[1] != 2.0) {
        return 1;
    }
    return matcal_linalg_multi_tu_matrix_value_b() == 3.0 ? 0 : 1;
}
