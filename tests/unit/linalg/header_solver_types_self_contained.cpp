#include "MatCal/Linalg/SolverTypes.hpp"

int matcal_linalg_header_solver_types_probe() {
    MatCal::Linalg::SolverOptions options;
    return options.valid() ? 0 : 1;
}
