# Pivoted LU User Guide

Dense partial-pivot LU is exposed as a reusable factorization.

```cpp
#include "MatCal/Linalg/DenseSolver.hpp"

MatCal::Linalg::DenseMatrix a{{0.0, 2.0}, {1.0, 3.0}};
auto factor = MatCal::Linalg::factorize_dense_partial_pivot(a);
if (factor.success()) {
    double det = factor.factorization.determinant();
}
```

The factorization records the row permutation and determinant sign. Its
mathematical contract is `P A = L U`; legacy no-pivot LU remains separate
because its old return type cannot represent `P`.

Reference:

- [Linalg Reference](../reference/linalg-v0.1.md)
- [Solver Result Contract](../reference/solver-result-contract.md)
