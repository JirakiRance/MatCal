#include "MatCal/Calculus/Calculus.hpp"

int matcal_calculus_header_self_contained() {
    auto result = MatCal::Calculus::forward_difference([](double x) { return x; }, 1.0, 1.0e-6);
    return result.success() ? 0 : 1;
}
