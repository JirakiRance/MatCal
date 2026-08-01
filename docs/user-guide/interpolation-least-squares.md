# Interpolation and Least Squares User Guide

`MatCal::Interpolation` owns interpolation data and rejects invalid nodes rather
than keeping dangling references.

```cpp
#include "MatCal/Interpolation/LinearInterpolator.hpp"

MatCal::Interpolation::LinearInterpolator line({{0.0, 1.0}, {2.0, 5.0}});
double y = line.evaluate(0.5);
```

`MatCal::LeastSquares` keeps the legacy polynomial least-squares capability on a
structured modern core. The first version still documents normal-equation
limitations and reports rank failures instead of silently returning a fake
success.

Reference:

- [Interpolation Reference](../reference/interpolation-v0.1.md)
- [LeastSquares Reference](../reference/least-squares-v0.1.md)
