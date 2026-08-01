# Polynomial and Roots User Guide

`MatCal::Polynomial` provides owning dense polynomial operations.

```cpp
#include "MatCal/Polynomial/Polynomial.hpp"

MatCal::Polynomial::Polynomial p{1.0, 0.0, 1.0};
double value = p.evaluate(3.0);
```

`MatCal::Roots` provides scalar root solvers with structured convergence
results.

```cpp
#include "MatCal/Roots/Roots.hpp"

auto root = MatCal::Roots::solve_bisection(
    [](double x) { return x * x - 2.0; }, 0.0, 2.0);
```

Reference:

- [Polynomial Reference](../reference/polynomial-v0.1.md)
- [Roots Reference](../reference/roots-v0.1.md)
