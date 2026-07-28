#include "MatCal/Nonlinear/Nonlinear.hpp"

int header_nonlinear_self_contained() {
    MatCal::Nonlinear::NonlinearOptions options;
    return options.valid() ? 0 : 1;
}
