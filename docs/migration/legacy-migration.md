# Legacy Migration

Legacy exists to keep historical MatCal users source-compatible. It is not a
signal that new code should keep using the old interfaces.

## What Remains

Installed legacy headers:

- `Basics.hpp`
- `Insert.hpp`
- `Iteration.hpp`
- `Matrix.hpp`
- `QinJiuShao.hpp`

Important namespaces and constructors are retained because PT and older
experiments may include them directly.

## Migration Direction

```text
legacy API -> compatibility facade -> modern core
```

New code should prefer:

- `MatCal::Linalg` for matrices and solvers;
- `MatCal::Polynomial` for polynomial operations;
- `MatCal::Roots` for scalar nonlinear equations;
- `MatCal::Interpolation` for interpolation;
- `MatCal::Calculus` and `MatCal::ODE` for numerical calculus and RK4;
- `MatCal::Nonlinear` and `MatCal::LeastSquares` for migrated higher-level algorithms.

See [Legacy Migration Matrix](legacy-migration-matrix.md) and
[Compatibility Policy](compatibility-policy.md).
