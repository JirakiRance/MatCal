# Numerical Contracts

M0 records legacy behavior and direction. M1 adds the first independent linalg contracts. M1.1 hardens linalg scale and finite-value behavior. M2 adds SPD skyline LDLT. This does not make every current algorithm numerically robust.

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

Prefer scale-aware tolerances:

```text
accept when residual <= max(abs_tol, rel_tol * scale)
```

The scale should be derived from matrix/vector norms and documented per solver.

## M1.1 Linalg Tolerance Contract

`MatCal::Linalg::SolverOptions` uses:

```text
max(absolute_tolerance, relative_tolerance * scale)
```

For pivot checks:

```text
pivot_factor * max(absolute_tolerance, relative_tolerance * scale)
```

Defaults:

- `absolute_tolerance = 0`
- `relative_tolerance = 1e-12`
- `pivot_factor = 1`
- `max_iterations = 1000`

The default absolute tolerance is zero so small-scale nonsingular systems are not rejected merely because their entries are below a fixed global threshold. The default relative tolerance is `1e-12`, chosen as a conservative double-precision reference value for the current dense baseline.

For pivot checks, the dense reference solver uses `scale = max(abs(A_ij))`, not `max(normInf(A), 1)`. The scale is not clamped to one. If `scale = 0`, the default tolerance is zero and an exactly zero pivot is reported as singular.

All numeric options must be finite; negative tolerances and negative `pivot_factor` are invalid; `max_iterations` must be nonzero. `pivot_factor = 0` is allowed and means exact-pivot mode. Tolerance arithmetic saturates to the largest finite double instead of producing Inf.

`Vector::norm2()` uses a scaled accumulation algorithm so values like `{1e200, 1e200}` do not overflow through an intermediate square sum.

`DenseMatrix` and `Vector` expose `all_finite()` checks. The dense reference solver rejects non-finite input with `SolverStatus::non_finite_input`.

`DenseMatrix::normInf()` can return Inf when a finite row sum overflows. The dense solver uses a separate safe maximum-absolute-entry scale so finite input does not become non-finite during scale computation without diagnosis.

Breakdown is distinct from singular:

- `singular`: a finite pivot is zero or below the pivot tolerance.
- `breakdown`: a finite-input computation produces NaN/Inf during factorization, back substitution, solution scaling, or residual evaluation.

M2.1 correction: the M1.1 implementation used `absolute_tolerance + relative_tolerance * scale`. M2.1 changes this to `max(absolute_tolerance, relative_tolerance * scale)` so callers can express `c * max(scale, 1)` by setting both tolerances to `c`.

Residual contract:

```text
abs_res = ||Ax-b||_inf
rel_res = abs_res / (matrix_scale * solution_scale + rhs_scale)
```

If the denominator is zero, `rel_res = 0` only when `abs_res = 0`. Success requires `abs_res <= max(absolute_tolerance, relative_tolerance * denominator)`.

## M2 Skyline LDLT Contract

`factorize_skyline_ldlt()` supports only real symmetric positive definite matrices. It performs an unpivoted LDLT factorization and rejects non-positive pivots with `SolverStatus::not_positive_definite`.

Status distinction:

- `not_positive_definite`: finite pivot is not positive above the pivot tolerance.
- `non_finite_input`: stored matrix values or RHS contain NaN/Inf.
- `breakdown`: finite input produces a non-finite intermediate during factorization or solve.
- `not_converged`: finite solve completes but residual exceeds the acceptance contract.

The skyline matrix scale is `max(abs(stored_value))`. Pivot tolerance reuses:

```text
pivot_factor * max(absolute_tolerance, relative_tolerance * matrix_scale)
```

Residual metrics use the same definition. The factorization owns its factors and can solve multiple right-hand sides without modifying the original matrix. M2.1 records `factorization_operation_count` and `solve_operation_count` separately, with `operation_count` retained as the total/reference counter.

## Solver Policy Direction

Automatic fallback, such as trying Gauss-Seidel after direct solve failure, should move out of low-level solvers into explicit policy functions. Low-level solver APIs should report why they failed.
