# MatCal::Calculus v0.1

M5 adds `MatCal::Calculus` as an independent header-only target for scalar finite differences and one-dimensional quadrature. It is not part of `MatCal::Linalg`, `MatCal::Roots`, or `MatCal::Polynomial`.

## Public Header

```cpp
#include <MatCal/Calculus/Calculus.hpp>
```

## Result Contract

`DerivativeResult` records:

- `value`;
- actual `step`;
- `function_evaluations`;
- machine-readable `CalculusDiagnostic`.

`IntegrationResult` records:

- `value`;
- `converged`;
- `estimated_error` when the algorithm has an estimate;
- iteration/refinement and function-evaluation counts;
- requested tolerance;
- optional Romberg table.

`success()` is true only for `CalculusStatus::success`. Invalid input, non-finite function values, and Romberg non-convergence are not reported as successful numeric results.

## Differentiation

`forward_difference` reuses the legacy `Derivative::dy_dx` formula:

```text
(f(x + h) - f(x)) / h
```

`central_difference` reuses the legacy centered formula:

```text
(f(x + h/2) - f(x - h/2)) / h
```

Both require a non-empty callable, finite `x`, finite positive `h`, finite shifted points, and finite function values. M5 does not implement automatic step-size selection.

## Integration

`integrate_instant` keeps the legacy left-rectangle rule but gives it a defined signed-interval contract. `a == b` returns exactly zero. For `a > b`, the integral is evaluated over the ordered interval and returned with a negative sign. The last rectangle is clipped to avoid stepping beyond the endpoint.

`integrate_newton_cotes` and `integrate_composite_newton_cotes` reuse the legacy closed Newton-Cotes coefficient table for orders 1 through 7. Reverse intervals use the same signed convention.

`integrate_romberg` reuses the legacy trapezoid-refinement and Richardson extrapolation loop shape. It now returns `not_converged` when the tolerance is not reached instead of returning the last estimated error as though it were an integral value.

## Legacy Delegation

The following legacy APIs now forward to this core:

- `Derivative::dy_dx`;
- `Derivative::dy_dx_center`;
- `NumericalIntegration::Instant`;
- `NumericalIntegration::NewtonCotes`;
- `NumericalIntegration::CompositeNewtonCotes`;
- `NumericalIntegration::Romberg`.

`Derivative::pF_px`, `Derivative::dF_dx`, least squares, and other `Basics.hpp` algorithms remain legacy-only.
