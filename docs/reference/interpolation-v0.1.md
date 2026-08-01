# MatCal Interpolation v0.1

Stage M4 introduced `MatCal::Interpolation` as an independent interpolation core. Stage M6 extends it to the remaining classic polynomial interpolation families. It provides general numerical interpolation only; it has no FEM, SFL, mesh, material, load, or result-recovery semantics.

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

## Polynomial Interpolation

M6 adds:

- `interpolate_lagrange`
- `interpolate_newton_divided`
- `interpolate_newton_finite`
- `interpolate_hermite`
- `DividedDifferenceResult`

The Lagrange core reuses the legacy Lagrange basis formula. Inputs must contain at least two finite points and no duplicate `x` values.

The Newton divided-difference core reuses the legacy lower-triangular divided-difference recurrence and exposes the table for compatibility diagnostics. Duplicate nodes are rejected before denominator division.

The Newton finite-difference core reuses the legacy forward-difference table and basis:

```text
sum_i delta^i y / (h^i i!) * product_{k=0}^{i-1}(x - (x0 + k h))
```

It requires finite `x0`, finite positive `h`, and at least two finite `y` values. Non-equidistant data should use divided differences instead.

The Hermite core reuses the legacy formula based on squared Lagrange basis functions:

```text
H_i(x) = (1 - 2 (x - x_i) l_i'(x_i)) l_i(x)^2
G_i(x) = (x - x_i) l_i(x)^2
```

Inputs must provide matching finite `xs`, `ys`, and derivative vectors. Nodes must be unique.

Legacy `LagrangeInsert`, `NewtonInsert_Quotient`, `NewtonInsert_Finite`, and `Hermite` now delegate to this core while preserving constructors, `calculate`, `reconstruct`, `getPoly`, and the Newton table getters.

## Still Not Covered

No barycentric interpolation, monotone spline, Akima spline, rational interpolation, or multidimensional interpolation API is provided in M6.
