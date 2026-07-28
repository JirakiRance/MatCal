# MatCal Nonlinear v0.1

Stage M6 adds `MatCal::Nonlinear` as an independent target for domain-neutral nonlinear systems. It does not depend on Legacy and does not contain SFL, FEM, mesh, material, load, or result semantics.

## Public Target

```cmake
target_link_libraries(app PRIVATE MatCal::Nonlinear)
```

`MatCal::Nonlinear` currently provides Newton iteration for square systems:

- `NonlinearStatus`
- `NonlinearReason`
- `NonlinearPhase`
- `NonlinearOptions`
- `NonlinearDiagnostic`
- `NonlinearMetrics`
- `NonlinearResult`
- `solve_newton_system`
- `solve_newton_system_finite_difference`

## Algorithm Source

The legacy `NewtonForEquations` implementation used the standard Newton equation:

```text
J(x_k) * delta = -F(x_k)
x_{k+1} = x_k + delta
```

M6 reuses that equation and the finite-difference Jacobian idea. The new core replaces the legacy Matrix solve and fallback policy with `MatCal::Linalg::solve_dense_partial_pivot` and structured diagnostics.

## Contract

Inputs must provide a finite initial guess and a residual of the same size. Analytic Jacobians must be square and match the residual dimension. The finite-difference path uses a caller-visible positive finite `finite_difference_step`.

Failures are structured:

- `invalid_input`: empty callable, invalid options, or empty state.
- `dimension_mismatch`: residual or Jacobian size mismatch.
- `non_finite_input`: non-finite initial state.
- `singular_jacobian`: dense linear solve reports a singular Jacobian.
- `breakdown`: residual, Jacobian, step, or update becomes non-finite.
- `not_converged`: maximum iterations are reached.

No partial solution is returned for numerical failure. Metrics record iterations, function evaluations, Jacobian evaluations, linear solves, final residual norm, final step norm, and requested tolerances.

## Legacy Facade

`MatCal::Algorithm::Iteration::NewtonForEquations::solve` keeps its public class, function name, argument order, and result struct. It now delegates to the finite-difference Newton core.

Old behavior could return the last iterate after a linear solve failure or divergence. M6 changes that failure path to `converged=false` with an empty `root`, preserving source compatibility while avoiding pseudo-success.
