# Linalg User Guide

Use `MatCal::Linalg` for modern matrix and solver work. It is independent of
`MatCal::Legacy`.

```cpp
#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/Vector.hpp"

MatCal::Linalg::DenseMatrix a{{4.0, 1.0}, {1.0, 3.0}};
MatCal::Linalg::Vector b{1.0, 2.0};
auto result = MatCal::Linalg::solve_dense_partial_pivot(a, b);
```

Always check `result.success()` before reading `result.solution`. Diagnostics
and metrics are machine-readable enough for callers to distinguish singular,
non-finite, non-converged, and breakdown cases without parsing stdout.

Reference:

- [Linalg Reference](../reference/linalg-v0.1.md)
- [Solver Result Contract](../reference/solver-result-contract.md)
- [Numerical Contracts](../reference/numerical-contracts.md)
