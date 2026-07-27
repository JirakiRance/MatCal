#include "MatCal/Linalg/Vector.hpp"

int matcal_linalg_header_vector_probe() {
    MatCal::Linalg::Vector v{1.0};
    return v.size() == 1 ? 0 : 1;
}
