# MatCal Interpolation v0.1

Stage M4 introduces `MatCal::Interpolation` as an independent interpolation core. It provides general numerical interpolation only; it has no FEM, SFL, mesh, material, load, or result-recovery semantics.

## Public Target

Use `MatCal::Interpolation`.

```cmake
target_link_libraries(app PRIVATE MatCal::Interpolation)
```

## LinearInterpolator

`LinearInterpolator` owns `xs` and `ys`. Inputs must be finite, same size, at least two points, and strictly increasing in `x`.

The segment formula is the legacy `LinearInsert` formula:

`y = y0 + (x - x0) / (x1 - x0) * (y1 - y0)`

Interval location uses binary search. Out-of-range behavior is explicit:

- `reject` is the new default.
- `clamp` returns the nearest endpoint value.
- `extrapolate` extends the nearest segment linearly.

Legacy `LinearInsert` delegates to this core with `extrapolate` to preserve old behavior.

## CubicSpline

`CubicSpline` implements the natural cubic spline currently used by legacy `CubicSpline`: endpoint second derivatives are zero, and interior second derivatives are solved from the tridiagonal system derived from adjacent slopes.

The new core reuses the legacy spline equation and evaluation formula but replaces the Legacy Matrix solve with an owned Thomas tridiagonal solve. It owns all nodes and second-derivative coefficients.

For two nodes, the spline degenerates to linear interpolation with zero endpoint second derivatives.

Legacy `CubicSpline` delegates construction and evaluation to this core with `extrapolate`, while preserving `calculate`, `getXs`, `getYs`, and `getM`.

## Not Migrated

`LagrangeInsert`, `NewtonInsert_Quotient`, `NewtonInsert_Finite`, and `Hermite` remain legacy-only in M4 except that they benefit indirectly from the migrated `QinJiuShao` polynomial core.
