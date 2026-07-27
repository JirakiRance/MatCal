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

- `absolute_tolerance = 1e-12`
- `relative_tolerance = 1e-12`
- `pivot_factor = 1.0`
- `max_iterations = 1000`

All numeric options must be finite. Negative tolerances are invalid. `max_iterations` must be nonzero.

Comparison form:

```text
absolute_tolerance + relative_tolerance * scale
```

Pivot comparison form:

```text
pivot_factor * (absolute_tolerance + relative_tolerance * scale)
```

For the M1 dense reference solver, `scale = max(normInf(A), 1)`.

## Diagnostics

`SolverDiagnostic` is deliberately domain-neutral:

- `status`
- `message`

It does not contain SFL diagnostics, source locations, AST nodes, materials, elements, loads, regions, JSON, or FEM result fields.

## Metrics

`SolverMetrics` currently contains:

- `iterations`
- `residual_norm`

For `solve_dense_partial_pivot`, `iterations` records elimination steps and `residual_norm` is the infinity norm of `Ax - b`.

## Result

`SolverResult` contains:

- `status`
- `solution`
- `metrics`
- `diagnostics`
- `method`
- `implementation`

Numerical failures should normally return a non-success `SolverResult`. Programmer errors and invalid object access, such as out-of-range indexing or impossible matrix allocation sizes, may still throw exceptions.
