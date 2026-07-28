#include "MatCal/Roots/Roots.hpp"

int matcal_roots_header_self_contained() {
    MatCal::Roots::RootOptions options;
    return options.valid() ? 0 : 1;
}
