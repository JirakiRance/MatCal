# Numerical Contracts

M0 records legacy behavior and direction. M1 adds the first independent linalg contracts. M1.1 hardens linalg scale and finite-value behavior. M2 adds SPD skyline LDLT. M4 adds scalar roots plus linear and natural cubic spline interpolation contracts. M5 adds scalar calculus and ODE contracts. M6 adds multivariable Newton systems, polynomial least squares, and remaining classic polynomial interpolation contracts. M7 adds dense stationary iterative linear solvers, dense power eigen solvers, and multivariable finite-difference derivative helpers. This does not make every current algorithm numerically robust.

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
- `NumericalIntegration::Instant` requires a non-empty callable, finite endpoints, and finite positive `eps`. It uses a left-rectangle rule and clips the final step to avoid integrating beyond the endpoint. M5 allows reverse intervals and returns the signed integral; `a == b` returns exactly zero.
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

## M4 Roots Contract

`MatCal::Roots` separates root status from the last iterate. `RootResult::success()` is true only when `RootStatus::success` and `converged` are both true.

Scalar root options use absolute and relative tolerances for iterate/step acceptance. Bisection accepts exact zero residual or interval-step convergence; it does not accept tiny absolute residual alone because function scaling is caller-defined and can be arbitrarily small.

Machine-readable diagnostics distinguish:

- invalid callable/options/input;
- interval or sign-change failure;
- derivative near zero;
- secant denominator near zero;
- non-finite function value or iterate;
- maximum iteration exhaustion.

Legacy `solveDetailed` functions map this contract back to old structs. Legacy throwing APIs still throw on invalid setup and numerical breakdown, and throw on non-convergence where they historically did.

## M4 Interpolation Contract

`MatCal::Interpolation::LinearInterpolator` and `CubicSpline` own their nodes. `xs` and `ys` must have matching size, at least two nodes, finite entries, and strictly increasing `x`.

`LinearInterpolator` uses binary interval lookup and the legacy segment formula. Extrapolation is explicit: reject, clamp, or extrapolate. The new default is reject; legacy `LinearInsert` uses extrapolate for source behavior compatibility.

`CubicSpline` is a natural cubic spline. Endpoint second derivatives are zero. The interior tridiagonal system and evaluation formula match the legacy `CubicSpline` implementation, but the new core uses an owned tridiagonal solve rather than Legacy Matrix dynamic polymorphism. For two nodes it degenerates to a line with zero second derivatives.

M6 polynomial interpolation functions own their input-derived coefficients through `MatCal::Polynomial::Polynomial`. Lagrange, divided-difference Newton, finite-difference Newton, and Hermite interpolation reject non-finite nodes, duplicate nodes, and invalid sizes before constructing a polynomial. Newton finite differences require a finite positive spacing `h`; non-equidistant data is outside that API contract.

## M5 Calculus Contract

`MatCal::Calculus` rejects empty callables, non-finite coordinates, non-finite step sizes, non-positive finite-difference steps, invalid Newton-Cotes orders, invalid segment counts, and non-finite function values. New APIs return structured `CalculusDiagnostic` data instead of printing.

Finite differences use the legacy formulas:

```text
forward: (f(x + h) - f(x)) / h
central: (f(x + h/2) - f(x - h/2)) / h
```

No automatic optimal step-size selection is claimed.

Integration contracts:

- zero-length interval returns exactly `0`;
- reverse interval returns the negative of the forward integral;
- left-rectangle `Instant` keeps the legacy left-endpoint rule;
- closed Newton-Cotes supports legacy orders 1 through 7;
- Romberg reports `not_converged` when the tolerance is not reached instead of returning an error estimate as a pseudo integral.

## M5 ODE Contract

`MatCal::ODE` accepts value-owned vector states and RHS callables of the form `f(t, y) -> vector`. It rejects empty states, non-finite time/state/step values, RHS size mismatches, and non-finite RHS outputs with structured diagnostics.

`rk4_step` uses the classic four-stage RK4 formula extracted from the legacy `ODE::RungeKutta_44` and PT-style `Integrate::RK4::step/step2` implementations. The new core allows negative `dt` for backward stepping by default; legacy table APIs still require positive `h`.

Trajectory rows are owned and use the compatibility shape:

```text
[t, y0, y1, ...]
```

## M6 Nonlinear Contract

`MatCal::Nonlinear` solves square vector systems with Newton iteration. It reuses the legacy Newton equation `J * delta = -F` but returns structured `NonlinearResult` values.

Options require finite non-negative absolute and relative tolerances, a finite positive finite-difference step, nonzero `max_iterations`, and valid linalg options.

Success requires a finite residual that satisfies:

```text
residual_norm <= max(abs_tol, rel_tol * max(initial_residual_norm, 1))
```

Step convergence is accepted only together with residual acceptance. Singular Jacobians are reported separately from non-finite breakdown. Numerical failures return no partial solution.

## M6 Least-Squares Contract

`MatCal::LeastSquares` keeps the legacy weighted normal-equation formula for polynomial fitting. This is a compatibility baseline, not a claim of high numerical stability.

Inputs require matching finite `x`, `y`, and weight vectors. Weights must be strictly positive in the M6 core. Selected degree lists must be non-empty and non-negative.

Rank-deficient normal equations return `LeastSquaresStatus::rank_deficient`. Overflow while building powers, normal equations, or the right-hand side returns `breakdown`. Successful results include coefficients, selected degrees, the fitted polynomial, normal matrix, right-hand side, sample count, term count, and infinity residual over the samples.

## M7 Stationary Iterative Solver Contract

`solve_jacobi`, `solve_gauss_seidel`, and `solve_sor` operate on `MatCal::Linalg::DenseMatrix` and `Vector`. They reuse the legacy stationary update formulas, with SOR using the corrected relaxed Gauss-Seidel update:

```text
x_new[i] = omega * gs_value[i] + (1 - omega) * x_old[i]
```

Inputs must be square, dimension-compatible, finite, and have finite diagonal entries above the diagonal tolerance:

```text
pivot_factor * max(absolute_tolerance, relative_tolerance * matrix_scale)
```

`omega` must be finite and satisfy `0 < omega < 2`. Success requires both a finite iterate and residual acceptance under the existing linalg residual contract. `not_converged` returns no partial solution in the core result. Metrics record iterations, operation count, matrix/RHS/solution scale, pivot tolerance, minimum diagonal magnitude, and absolute/relative residuals.

M7.1 clarification: residual evaluation is part of the numerical process. If a finite matrix, RHS, and initial guess produce a non-finite initial residual, the solver returns `SolverStatus::breakdown` immediately and does not continue iterating.

## M7 Dense Power Eigen Contract

`dominant_eigenpair` and `inverse_power_eigenpair` are dense reference algorithms in `MatCal::Linalg`. They are not sparse eigensolvers and do not claim robust behavior for clustered, repeated, or defective spectra.

Power iteration normalizes by the largest-magnitude vector component, computes the eigenvalue by Rayleigh quotient, and accepts convergence only when the eigen residual is finite and satisfies:

```text
||A v - lambda v||_inf <= max(abs_tol, rel_tol * max(matrix_scale * vector_scale, |lambda| * vector_scale, 1))
```

Inverse power applies shifted solves with `solve_dense_partial_pivot(A - shift I, rhs)`. Singular shifted systems return `EigenStatus::singular_shift`; non-finite products, quotients, or solve intermediates return `breakdown`. Failures do not return a pseudo-success eigenpair.

## M7 Multivariable Derivative Contract

`MatCal::Calculus::partial_difference` and `gradient` reuse the legacy forward finite-difference formula:

```text
partial_i F(x) = (F(x + h e_i) - F(x)) / h
```

The helpers require a non-empty callable, finite owned input state, an in-range coordinate, and finite positive step. They do not mutate caller state. They return structured derivative results with evaluation counts and reject non-finite function values instead of propagating them silently.
