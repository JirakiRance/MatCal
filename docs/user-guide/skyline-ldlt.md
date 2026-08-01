# Skyline LDLT User Guide

`SymmetricSkylineMatrix` stores a real symmetric matrix by row profile. The
current factorization is no-pivot LDLT for symmetric positive definite systems.

```cpp
#include "MatCal/Linalg/SkylineLdlt.hpp"
#include "MatCal/Linalg/SymmetricSkylineMatrix.hpp"

MatCal::Linalg::SymmetricSkylineMatrix a({0, 0});
a.set(0, 0, 4.0);
a.set(1, 0, 1.0);
a.set(1, 1, 3.0);

auto factor = MatCal::Linalg::factorize_skyline_ldlt(a);
```

The factorization owns its data and can solve multiple right-hand sides without
re-factorizing. MatCal only sees numeric matrix/profile data; FEM assembly,
constraints, loads, and result recovery stay in downstream projects such as SFL.

Reference:

- [Skyline Storage and SPD LDLT](../reference/skyline-ldlt.md)
- [SFL Integration Boundary](../reference/sfl-integration-boundary.md)
