# Known Limitations and Confirmed Bugs

This list is intentionally direct. Do not treat confirmed bugs as the standard for future new APIs.

## Confirmed Bugs / High-Risk Behavior

- `NumericalIntegration::Instant` assigns `_func = nullptr` instead of checking it. Risk: any valid function is cleared, leading to runtime failure when called. Compatibility handling: add a safe replacement, then keep legacy wrapper with a bug note or corrected behavior at a version boundary.
- `SOR` takes `int omega`, preventing normal relaxation factors such as `1.25`, and uses `dense_A->get(i,k)` where the old iterate should be used. Risk: wrong results except cases like `omega == 1` where the bug is masked. Compatibility handling: add `SOR(..., double omega, ...)` or a new solver API; keep legacy int overload deprecated.
- `LU_Decompose` has no pivoting. Risk: fails on invertible matrices with zero/small leading pivots. Compatibility handling: document as no-pivot LU; add pivoted LU under a new name/result type.
- `solve_columnElimination` lacks complete dimension and pivot validation. Risk: singular systems can produce division by zero or unclear failures. Compatibility handling: add explicit status-based direct solver.
- `determinant` ignores row-swap sign from column elimination. Risk: determinant sign can be wrong after pivot swaps. Compatibility handling: fix in new implementation; preserve old behavior only through legacy if needed.
- `Matrix(AbstractMatrix&)` and some triangular copy constructors do not explicitly initialize base dimensions. Risk: copied objects can have mismatched metadata and storage. Compatibility handling: add tests before fixing; fix as source-compatible bug repair if no external code depends on broken metadata.
- `TridiagonalMatrix(0)` routes through `resize(0,0)` and resizes `n - 1`. Risk: negative-to-size_t conversion may attempt huge allocation. Compatibility handling: fix zero-size handling with tests.
- `OrthogonalPolynomials::Legendre` uses integer division in recurrence coefficients. Risk: wrong polynomial coefficients for higher orders. Compatibility handling: add numeric tests and fix in M1 or M2.
- `QinJiuShao::toFunction()` captures `this`. Risk: dangling callable after polynomial destruction. Compatibility handling: add `toOwningFunction()` or change new API to capture a value; deprecate borrowed callable.

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
