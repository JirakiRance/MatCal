# MatCal Roots v0.1

Stage M4 introduces `MatCal::Roots` as an independent scalar nonlinear equation core. It does not depend on Legacy, Linalg, SFL, FEM, or application diagnostics.

## Public Target

Use `MatCal::Roots`.

```cmake
target_link_libraries(app PRIVATE MatCal::Roots)
```

## Migrated Algorithms

- `solve_bisection`
- `solve_picard`
- `solve_picard_aitken`
- `solve_newton`
- `solve_downhill_newton`
- `solve_secant_two_point`
- `solve_secant_one_point`

The formulas and iteration skeletons are taken from `src/Iteration.hpp`: bisection interval halving, Picard fixed-point iteration, Aitken acceleration, Newton step, downhill Newton half-step shrink, and secant update formulas.

## Result Contract

`RootResult` returns the candidate value only with a structured status:

- `success` means convergence was accepted.
- `not_converged` means no root is claimed.
- derivative-zero, denominator-zero, bracket failure, invalid input, and non-finite intermediate states are machine-readable.

Bisection accepts convergence by exact zero function value or interval-step tolerance. It does not treat a tiny absolute residual alone as proof of a root, because scaled functions such as `1e-20 * (x - c)` would otherwise accept arbitrary endpoints.

## Legacy Delegation

The following legacy APIs now delegate to `MatCal::Roots` and map the result back to old result structs or exceptions:

- `MatCal::Algorithm::Iteration::Bisection`
- `Picard`
- `Picard::solve_Aitken`
- `Newton`
- `Newton::solve_downhill`
- `Secant`

`NewtonForEquations` is not migrated in M4.
