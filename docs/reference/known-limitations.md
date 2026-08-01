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

- `LU_Decompose` has no pivoting. Risk: fails on invertible matrices with zero/small leading pivots. Compatibility handling: document as no-pivot LU; use `MatCal::Linalg::factorize_dense_partial_pivot` when a reusable pivoted factorization is required.
- `solve_columnElimination` now delegates through the structured dense partial-pivot LU factorization, but its legacy facade still maps numerical failure to exceptions rather than returning `SolverResult`.
- `solve_Linear_System` automatically tries multiple algorithms. Risk: policy and diagnostics are mixed into a low-level helper. Compatibility handling: preserve legacy behavior, then add explicit status-based solver selection.
- `Matrix(AbstractMatrix&)` still takes a non-const reference even though it copies. Risk: const callers cannot use it directly. Compatibility handling: add const-safe overloads in M1 without removing the old constructor.
- Several dynamic casts assume dense conversion success. Risk: future matrix types could expose unchecked null paths. Compatibility handling: centralize checked conversions in M1 internals.
- Some numerical entry points still accept NaN/Inf without explicit contracts.

## Design Limitations

- Header-only definitions previously had ODR/multiple-definition risk; M0 marks known non-template definitions `inline`.
- Many APIs use non-const references for inputs that are not intentionally mutated.
- Many dynamic casts assume `toNormalMatrix()` returns `Matrix`; several are checked, but not all.
- Core code prints to stdout in legacy `show()` helpers. M6 removed the `NewtonForEquations` validation-path stdout from the migrated solve path.
- Sparse storage is triplet/vector based with linear lookup; duplicate triplets are not merged by the constructor.
- Some remaining legacy-only paths still use fixed absolute tolerances.
- No thread-safety contract exists. Independent objects are generally usable independently; shared mutable objects are not synchronized.
- No ABI stability contract exists.
- No Windows-only `windows.h` dependency was found in the core headers during M0.

## M1 Linalg Limitations

- `MatCal::Linalg` is a 0.x development API, not an ABI-stable API.
- `DenseMatrix` is contiguous row-major storage for small dense work and reference tests. It is not a replacement for future large FEM sparse storage.
- `solve_dense_partial_pivot` is a one-shot facade over reusable partial-pivot LU. It is not optimized for large systems. Use `factorize_dense_partial_pivot` plus `PivotedLuFactorization::solve(DenseMatrix)` for multiple right-hand sides.
- `not_positive_definite` is used by skyline LDLT and remains available for future SPD algorithms. The dense reference solver can use `not_converged` if a finite residual exceeds the acceptance contract.
- M1 does not bridge legacy APIs to `MatCal::Linalg`; the targets intentionally remain parallel.

## M1.1 Fixed Linalg Contract Issue

- M1 used `max(normInf(A), 1)` and `absolute_tolerance = 1e-12` for pivots. This could reject `1e-20 * I` even though it is nonsingular and condition number 1.
- M1.1 changes the default absolute tolerance to `0`, keeps relative tolerance at `1e-12`, and uses an unclamped maximum-absolute-entry matrix scale for pivot checks.
- Finite intermediate overflow in the dense solver is now reported as `breakdown`, not `singular`.

## M2 Linalg Limitations

- `SymmetricSkylineMatrix` stores a symmetric numerical profile only. It does not perform FEM assembly or constraint handling.
- `SkylineLdltFactorization` is SPD-only and unpivoted. It rejects singular, semidefinite, and indefinite matrices with `not_positive_definite`.
- General indefinite symmetric factorization, Bunch-Kaufman pivoting, CSR, and iterative methods remain future work.

## M2.1 Fixed Linalg Contract Issue

- M2 skyline back substitution scanned all later rows and filtered by profile. M2.1 stores column adjacency in the factorization so back substitution visits only profile rows containing the active column.
- M2.1 changes linalg tolerance composition from `absolute_tolerance + relative_tolerance * scale` to `max(absolute_tolerance, relative_tolerance * scale)`, matching the integration need to express `c * max(scale, 1)`.

## M3 Polynomial Migration Notes

- `QinJiuShao` no longer owns independent implementations of Horner evaluation, arithmetic, derivative, integral, definite integral, or owning callable creation. Those paths delegate to `MatCal::Polynomial::Polynomial`.
- The legacy `QinJiuShaoNode` negative-degree constructor no longer prints to stdout before throwing.
- `MatCal::Polynomial::Polynomial` uses dense coefficient storage. M4 adds an explicit dense-degree limit so very high sparse terms fail before large allocation attempts. A future sparse polynomial representation should be added only after measured need.

## M4 Roots and Interpolation Migration Notes

- Legacy scalar root classes now delegate to `MatCal::Roots`.
- Bisection no longer accepts a tiny absolute residual alone as endpoint success. This fixes scaled-function false positives.
- Legacy `LinearInsert` and `CubicSpline` now delegate to `MatCal::Interpolation`.
- M6 migrates `LagrangeInsert`, `NewtonInsert_Quotient`, `NewtonInsert_Finite`, and `Hermite`.

## M5 Calculus and ODE Migration Notes

- Legacy scalar `Derivative::dy_dx` and `Derivative::dy_dx_center` now delegate to `MatCal::Calculus`; M7 extends this to `pF_px` and `dF_dx`.
- Legacy `NumericalIntegration::Instant`, `NewtonCotes`, `CompositeNewtonCotes`, and callable `Romberg` now delegate to `MatCal::Calculus`.
- `Instant` now supports reverse intervals with signed results. Code that depended on reverse intervals throwing should check intervals before calling.
- `Romberg` now reports non-convergence through the legacy throwing path instead of returning an error estimate as a successful integral value.
- Legacy `ODE::SimpleEuler`, `ODE::Euler`, `ODE::RungeKutta_44`, and PT-style `Integrate::RK4::step/step2` now delegate to `MatCal::ODE`.
- Legacy ODE table APIs still require positive `h`; the new `MatCal::ODE` core allows negative `dt` by option for backward stepping.
- M7 migrates `pF_px` and `dF_dx`.

## M6 Nonlinear, Least-Squares, and Interpolation Migration Notes

- Legacy `NewtonForEquations::solve` now delegates to `MatCal::Nonlinear`. It no longer returns the last finite iterate as a pseudo root on singular-Jacobian or numerical failure; the legacy result has `converged=false` and an empty `root`.
- `MatCal::Nonlinear` currently provides undamped Newton only. No line search, trust region, sparse Jacobian, or overdetermined system solve is provided.
- Legacy `Least_Square::solve` overloads now delegate to `MatCal::LeastSquares`. The core keeps normal equations for compatibility, so ill-conditioned fits remain a numerical limitation. QR/SVD are future work.
- M6 least-squares weights are strictly positive. Code that wants zero-weight masking should filter samples or use a future explicitly documented masking API.
- Legacy `LagrangeInsert`, `NewtonInsert_Quotient`, `NewtonInsert_Finite`, and `Hermite` now delegate to `MatCal::Interpolation`.
- Duplicate or non-finite interpolation nodes are rejected by the new core. Code that accidentally depended on duplicate-node division by zero must clean inputs before calling.

## M7 Matrix, Iterative, Eigen, and Derivative Migration Notes

- Legacy matrix-to-linalg adapters deep-copy data and reject non-finite values. They do not provide zero-copy views.
- Legacy `solve_columnElimination` now delegates to `MatCal::Linalg::factorize_dense_partial_pivot` and solves all right-hand-side columns with a single factorization.
- Legacy `Jacobi`, `Gauss_Seidel`, and `SOR` now delegate to `MatCal::Linalg` stationary solvers. Core non-convergence returns `SolverStatus::not_converged` without a partial solution; the legacy facade maps this to `converged=false`.
- `MatCal::Linalg` stationary solvers are dense reference algorithms. They do not include sparse matrix support, preconditioners, or convergence guarantees for non-diagonally-dominant systems.
- M7.1 fixes a stationary-solver edge case where an overflowed initial residual could continue into the iteration loop instead of immediately returning `breakdown`.
- Legacy `PowerMethod` and `PowerMethod_reverse` now delegate to `MatCal::Linalg` eigen solvers. Shifted inverse power uses the dense partial-pivot solver instead of legacy no-pivot LU.
- The M7 eigen solvers are basic dense power iterations. Repeated, clustered, or difficult eigenvalues may converge slowly or report `not_converged`.
- Legacy `LU_Decompose` remains no-pivot LU and is not silently replaced with pivoted LU.
- `Derivative::pF_px` and `Derivative::dF_dx` now delegate to `MatCal::Calculus` partial derivative and gradient helpers.

## M8 Dense Pivoted LU Migration Notes

- `MatCal::Linalg::PivotedLuFactorization` owns compact `LU`, row permutation, permutation sign, the original matrix copy used for residual checks, and factorization metrics.
- `solve_dense_partial_pivot` now delegates to `factorize_dense_partial_pivot` and then `PivotedLuFactorization::solve`.
- Dense multi-RHS solves are atomic and reuse one factorization.
- Legacy `determinant` now delegates to pivoted LU and computes `permutation_sign * product(diag(U))`.
- Legacy `LU_Decompose` remains the old no-pivot compatibility API because its result type cannot represent `P` in `P A = L U`.

## M8.1 Release Candidate Limitations

- `v0.3.0-alpha.1` is a source-package release candidate, not an ABI-stable release.
- The CMake package version is numeric `0.3.0`; the annotated Git tag should carry the prerelease name `v0.3.0-alpha.1` after remote CI is green.
- `MatCal::Linalg` still provides dense reference algorithms and SPD skyline LDLT only. No CSR, QR, SVD, Bunch-Kaufman, or FEM-specific solver is included.
- The SFL-compatible consumer is a generic API compatibility check. It is not SFL integration and does not validate FEM assembly, constraints, material behavior, result recovery, or SFL diagnostics.
- CI workflows are added for GitHub Linux GCC and Windows MSVC, but they should not be reported as passed until they actually run on the remote.
- MSVC add-subdirectory consumers that include legacy root headers may still report C4819 code-page warnings from historical non-ASCII comments in `src/Insert.hpp`, `src/QinJiuShao.hpp`, and `src/Matrix.hpp`. M8.1 does not perform a broad source encoding conversion.
