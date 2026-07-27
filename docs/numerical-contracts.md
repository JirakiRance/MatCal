# Numerical Contracts

M0 records current behavior and the desired direction. It does not make every current algorithm numerically robust.

## Current Legacy Contracts

- Dimensions are generally checked for matrix construction, indexing, and many operations.
- Most solvers use fixed absolute tolerances such as `1e-12` or `1e-6`.
- NaN and Inf are not rejected at boundaries; tests characterize propagation through simple matrix operations.
- Gaussian elimination copies inputs and uses partial row pivoting, but lacks complete dimension and zero-pivot checks.
- LU is Doolittle without pivoting and is valid only when all pivots are safely nonzero. A zero pivot can mean "pivoting required", not necessarily singular.
- Thomas/chase solver is valid only for suitable nonsingular tridiagonal systems without zero pivots.
- Iterative solvers report `converged`, `iterations`, and `error`, but convergence criteria are fixed absolute deltas.

## M0.1 Clarified Contracts

- `TridiagonalMatrix(0)` is allowed. Its lower, diagonal, and upper arrays are empty; solving against a `0 x k` right-hand side returns a `0 x k` matrix.
- `NumericalIntegration::Instant` requires a non-empty callable, finite endpoints, finite positive `eps`, `a <= b`, and `eps <= b - a`. It uses a left-rectangle rule and clips the final step to avoid integrating beyond `b`.
- `OrthogonalPolynomials::Legendre(n)` uses the standard recurrence and is regression-tested against analytic `P0`, `P1`, `P2`, and `P3` values.
- `determinant` copies its input, uses the existing row-pivoting elimination path, applies row-swap parity, and keeps the legacy singular path as an exception.
- `QinJiuShao::toFunction()` captures coefficients by value. The callable observes the polynomial state at creation time and may outlive the original object.
- `SOR` validates finite `omega` with `0 < omega < 2`. Convergence tests now check both reported iteration state and `Ax - b` residual.

## Desired New Contracts

Future APIs should define:

- Input shape and storage requirements.
- Whether input is mutated.
- Whether NaN/Inf is accepted, rejected, or propagated.
- Scale-aware absolute/relative tolerance.
- Residual norm used for convergence.
- Maximum iterations and stopping reason.
- Singular, ill-conditioned, non-converged, and invalid-option statuses.

## Tolerance Direction

Prefer combined tolerances:

```text
accept when residual <= abs_tol + rel_tol * scale
```

The scale should be derived from matrix/vector norms and documented per solver.

## Solver Policy Direction

Automatic fallback, such as trying Gauss-Seidel after direct solve failure, should move out of low-level solvers into explicit policy functions. Low-level solver APIs should report why they failed.
