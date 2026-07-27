# Numerical Contracts

M0 records current behavior and the desired direction. It does not make every current algorithm numerically robust.

## Current Legacy Contracts

- Dimensions are generally checked for matrix construction, indexing, and many operations.
- Most solvers use fixed absolute tolerances such as `1e-12` or `1e-6`.
- NaN and Inf are not rejected at boundaries; tests characterize propagation through simple matrix operations.
- Gaussian elimination copies inputs and uses partial row pivoting, but lacks complete dimension and zero-pivot checks.
- LU is Doolittle without pivoting and is valid only when all pivots are safely nonzero.
- Thomas/chase solver is valid only for suitable nonsingular tridiagonal systems without zero pivots.
- Iterative solvers report `converged`, `iterations`, and `error`, but convergence criteria are fixed absolute deltas.

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
