# Solver Result Contract

M1 introduces structured solver reporting in `include/MatCal/Linalg/SolverTypes.hpp`.

## Status

`SolverStatus` currently distinguishes:

- `success`
- `invalid_input`
- `dimension_mismatch`
- `non_finite_input`
- `singular`
- `not_positive_definite`
- `breakdown`
- `not_converged`

Not every solver uses every status. Unused statuses are reserved so direct and iterative solvers can share one vocabulary without adding SFL or FEM diagnostics.

## Options

`SolverOptions` contains:

- `absolute_tolerance = 0`
- `relative_tolerance = 1e-12`
- `pivot_factor = 1.0`
- `max_iterations = 1000`

All numeric options must be finite. Negative tolerances are invalid. `max_iterations` must be nonzero. `pivot_factor = 0` is valid and means exact-pivot mode: only an exactly zero pivot fails the pivot tolerance test.

The default absolute tolerance is zero so it does not dominate small-scale but well-conditioned systems such as `1e-20 * I`. The default relative tolerance is `1e-12`, a conservative multiple of double precision machine epsilon for the current reference solver. Future solvers may expose algorithm-specific options, but they should keep the same dimensional meaning.

Comparison form:

```text
absolute_tolerance + relative_tolerance * scale
```

Pivot comparison form:

```text
pivot_factor * (absolute_tolerance + relative_tolerance * scale)
```

M1.1 correction: for pivot checks the dense reference solver uses:

```text
matrix_scale = max(abs(A_ij))
```

`matrix_scale` is not clamped to `1`. If `matrix_scale = 0`, the default pivot tolerance is zero and an exactly zero pivot is reported as `singular`. Tolerance arithmetic uses controlled saturation to `double::max()` rather than producing Inf.

## Diagnostics

`SolverDiagnostic` is deliberately domain-neutral:

- `status`
- `code`
- `reason`
- `phase`
- `row`
- `column`
- `value`
- `scale`
- `tolerance`
- `message`

It does not contain SFL diagnostics, source locations, AST nodes, materials, elements, loads, regions, JSON, or FEM result fields.

`code`, `reason`, and `phase` are stable machine-readable fields. `message` is for humans and should not be parsed by consumers.

Current solver examples:

- `pivot_too_small` with status `singular`.
- `non_positive_pivot` with status `not_positive_definite`.
- `non_finite_intermediate` with status `breakdown`.
- `residual_too_large` with status `not_converged`.
- `invalid_options`, `matrix_not_square`, `rhs_size_mismatch`, and `non_finite_input` during validation.

## Metrics

`SolverMetrics` currently contains:

- `iterations`
- `operation_count`
- `residual_norm`
- `absolute_residual_norm`
- `relative_residual_norm`
- `residual_acceptance_tolerance`
- `matrix_scale`
- `rhs_scale`
- `solution_scale`
- `pivot_tolerance_used`
- `minimum_abs_pivot`

For `solve_dense_partial_pivot`, `iterations` records elimination steps and `residual_norm` is retained as an alias for `absolute_residual_norm`. For skyline LDLT, `iterations` records factorization rows or solve size depending on phase, and `operation_count` is a reference work counter rather than a hardware-level flop count.

Residual contract:

```text
absolute_residual_norm = ||Ax - b||_inf
relative_residual_norm = ||Ax - b||_inf / (matrix_scale * solution_scale + rhs_scale)
```

If the denominator is zero, relative residual is `0` when the absolute residual is zero. A nonzero residual with zero denominator is a numerical breakdown.

Success requires:

```text
absolute_residual_norm <= absolute_tolerance + relative_tolerance * (matrix_scale * solution_scale + rhs_scale)
```

Otherwise the result is non-success; the dense reference solver uses `not_converged` for finite residuals above tolerance and `breakdown` for non-finite residual evaluation.

## Result

`SolverResult` contains:

- `status`
- `solution`
- `metrics`
- `diagnostics`
- `method`
- `implementation`

Numerical failures should normally return a non-success `SolverResult`. Programmer errors and invalid object access, such as out-of-range indexing or impossible matrix allocation sizes, may still throw exceptions.

On non-success numerical results, dense and skyline solvers do not expose partial solutions.
