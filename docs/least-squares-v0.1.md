# MatCal LeastSquares v0.1

Stage M6 adds `MatCal::LeastSquares` as an independent target for polynomial least-squares fitting. It is a general numerical module and does not contain FEM or SFL semantics.

## Public Target

```cmake
target_link_libraries(app PRIVATE MatCal::LeastSquares)
```

Public API:

- `LeastSquaresStatus`
- `LeastSquaresReason`
- `LeastSquaresDiagnostic`
- `LeastSquaresMetrics`
- `LeastSquaresResult`
- `fit_polynomial`
- `fit_polynomial_degree`
- `fit_polynomial_selected`

## Algorithm Source

The M6 core reuses the legacy `Least_Square` construction of weighted normal equations over polynomial basis terms:

```text
A_ij = sum_k w_k * x_k^(degree_i + degree_j)
b_i  = sum_k w_k * y_k * x_k^degree_i
```

The selected-term overload preserves the old selected polynomial basis behavior. The output polynomial is now a `MatCal::Polynomial::Polynomial`, while legacy overloads convert it back to `QinJiuShao`.

## Contract

Inputs must have matching non-empty `x`, `y`, and `weights`. All values must be finite. M6 freezes weights as strictly positive for this first structured API. Selected degrees must be non-negative and non-empty.

The first implementation intentionally keeps normal equations. This is compatible with the old algorithm but remains less stable than QR or SVD for ill-conditioned problems. Rank deficiency is detected through the linalg solve and reported as `LeastSquaresStatus::rank_deficient`; it is not returned as a successful fit.

## Legacy Facade

All legacy `Least_Square::solve` overloads keep their public signatures and now delegate to `MatCal::LeastSquares`. The legacy degree rule `degree > 0` remains in the facade, while the new core can fit a degree-0 constant.

Old failure paths could return a placeholder polynomial. M6 keeps the old result struct but returns an empty coefficient vector and diagnostic message for structured failures.
