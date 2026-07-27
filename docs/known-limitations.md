# Known Limitations and Confirmed Bugs

This list is intentionally direct. Do not treat confirmed bugs as the standard for future new APIs.

## M0.1 Fixed High-Risk Bugs

- `NumericalIntegration::Instant` now checks `_func` instead of assigning it to `nullptr`. Empty callables and non-finite inputs throw `std::invalid_argument`.
- `SOR` now has a `double omega` overload, validates finite `0 < omega < 2`, and uses the correct old iterate in the relaxation update. The legacy `int` overload remains as a wrapper.
- `determinant` now applies row-swap parity and keeps the input matrix unchanged.
- `Matrix` triangular copy/assignment behavior now preserves base dimensions and storage metadata.
- `TridiagonalMatrix(0)` now creates empty bands instead of attempting a huge `size_t` allocation.
- `OrthogonalPolynomials::Legendre` now uses floating-point recurrence coefficients.
- `QinJiuShao::toFunction()` now captures coefficient state by value, so returned callables may outlive the source polynomial.

## Remaining Confirmed Bugs / High-Risk Behavior

- `LU_Decompose` has no pivoting. Risk: fails on invertible matrices with zero/small leading pivots. Compatibility handling: document as no-pivot LU; add pivoted LU under a new name/result type.
- `solve_columnElimination` lacks complete dimension and pivot validation. Risk: singular systems can produce division by zero or unclear failures. Compatibility handling: add explicit status-based direct solver.
- `solve_Linear_System` automatically tries multiple algorithms. Risk: policy and diagnostics are mixed into a low-level helper. Compatibility handling: preserve legacy behavior, then add explicit status-based solver selection.
- `Matrix(AbstractMatrix&)` still takes a non-const reference even though it copies. Risk: const callers cannot use it directly. Compatibility handling: add const-safe overloads in M1 without removing the old constructor.
- Several dynamic casts assume dense conversion success. Risk: future matrix types could expose unchecked null paths. Compatibility handling: centralize checked conversions in M1 internals.
- Some numerical entry points still accept NaN/Inf without explicit contracts.

## Design Limitations

- Header-only definitions previously had ODR/multiple-definition risk; M0 marks known non-template definitions `inline`.
- Many APIs use non-const references for inputs that are not intentionally mutated.
- Many dynamic casts assume `toNormalMatrix()` returns `Matrix`; several are checked, but not all.
- Core code prints to stdout in `show()`, negative `QinJiuShaoNode`, and `NewtonForEquations` validation paths.
- Sparse storage is triplet/vector based with linear lookup; duplicate triplets are not merged by the constructor.
- Fixed absolute tolerances dominate numerical decisions.
- No thread-safety contract exists. Independent objects are generally usable independently; shared mutable objects are not synchronized.
- No ABI stability contract exists.
- No Windows-only `windows.h` dependency was found in the core headers during M0.

## M1 Linalg Limitations

- `MatCal::Linalg` is a 0.x development API, not an ABI-stable API.
- `DenseMatrix` is contiguous row-major storage for small dense work and reference tests. It is not a replacement for future large FEM sparse storage.
- `solve_dense_partial_pivot` is a reference Gaussian-elimination solver. It is not optimized for large systems and does not provide multiple right-hand sides.
- `not_positive_definite` and `not_converged` are status vocabulary entries for future algorithms; the M1 dense direct solver does not use them.
- Skyline, LDLT, CSR, and iterative solver infrastructure remain future work.
- M1 does not bridge legacy APIs to `MatCal::Linalg`; the targets intentionally remain parallel.
