# MatCal Code Audit

Audit date: 2026-07-27.

## Repository State

- Branch: `main`.
- HEAD: `ce295a7a39cbbd21c1a39dde4daedf9180aa2a0d` (`添加了readme`).
- Remotes:
  - `github`: `https://github.com/JirakiRance/MatCal.git`
  - `origin`: `https://gitee.com/jirakirance/mat-cal.git`
- Initial status: clean relative to `github/main`. Git emitted permission warnings for `C:\Users\DELL/.config/git/ignore`; repository state itself was clean.
- Tracked files at audit start: `.gitignore`, `README.md`, `src/Basics.hpp`, `src/Insert.hpp`, `src/Iteration.hpp`, `src/Matrix.hpp`, `src/QinJiuShao.hpp`.
- Untracked files at audit start: none reported by `git status --short --branch`.
- EOL baseline: `git ls-files --eol` reported `i/lf w/crlf` for every tracked file. This means the index stores LF while the working tree uses CRLF. No pure line-ending changes should be committed.

## Recent History

Recent commits show MatCal is an existing handwritten numerical library, not a greenfield project:

- `ce295a7` add README.
- `fbfc8c2` ignore build/editor outputs and stop tracking them.
- `bc2cdb9` end-of-course snapshot and bug fixes.
- `d1b9fae` Matrix parent-pointer constructor; orthogonal polynomials, numerical integration, ODE.
- `1f1c1bb` PT project updated spline interpolation and linear interpolation; add integration, RK4, Euler integration.
- `a2877c7` move least squares into `Least_Square::solve`.
- `a734a50` implement `TridiagonalMatrix` and chase/Thomas solver.
- Earlier commits add nonlinear equations, interpolation namespaces, SOR, Jacobi, Gauss-Seidel, LU and QinJiuShao.

## README Findings

`README.md` is an API manual for the current legacy library. It documents:

- `MatCal::Utils`: `QinJiuShaoNode`, `QinJiuShao`, `AbstractMatrix`, `Matrix`, `SparseMatrix`, upper/lower triangular matrices, `TridiagonalMatrix`.
- `MatCal::Algorithm::Matrix`: elementary row/column transformations, norms, Gaussian elimination, LU, Jacobi, Gauss-Seidel, SOR, power methods.
- `MatCal::Algorithm::Basics`: derivatives, least squares, orthogonal polynomials, numerical integration, ODE.
- `MatCal::Algorithm::Insert`: Lagrange, Newton quotient, Newton finite-difference, Hermite, cubic spline, linear interpolation.
- `MatCal::Algorithm::Iteration`: scalar and vector nonlinear solvers.

The README does not explicitly describe SFL. The Git history and source comments explicitly mention PT, especially interpolation and `MatCal::Algorithm::Basics::Integrate::RK4`.

## Build Baseline Added

This M0 pass added standard CMake support:

- `MatCal::Legacy`: current header-only legacy API target.
- C++20 compile features.
- CTest tests in `tests/`.
- GCC/Clang/MSVC warning options for test targets.
- Optional `MATCAL_ENABLE_SANITIZERS` for GCC/Clang-family builds.

No public header was removed or renamed.

## Actual Build/Test Result

Observed platform:

- OS shell: Windows PowerShell in `E:\Ducuments\codeForVScode\MatCal`.
- CMake: 4.4.0.
- Compiler used by CMake: GNU C++ 15.1.0, `E:/minGW/mingw64/bin/c++.exe`.
- Clang/MSVC: not found in PATH during this audit.

Commands run:

```powershell
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Debug
cmake --build build --config Debug
ctest --test-dir build --output-on-failure
```

CTest result: 2/2 tests passed.

Warnings remain. They are documented rather than batch-fixed in M0 because many are tied to legacy copy rules, signedness, and numerical contracts.

## Engineering Changes Made

- Added `CMakeLists.txt`.
- Added characterization tests.
- Added multi-translation-unit link test.
- Marked header-defined non-template functions/static data as `inline` where required to make the legacy header-only library link safely from multiple translation units.
- Added missing standard includes needed by the existing header definitions.

No algorithmic rewrite was performed.

## M0.1 Safety Baseline

M0.1 keeps the legacy API surface intact and applies targeted source-compatible safety repairs only after adding regression coverage. It does not add SFL integration, FEM concepts, CSR, Skyline, LDLT, or a stable `MatCal::Linalg` target.

New regression coverage:

- Triangular matrix copy, assignment, transpose, solve, and metadata.
- Zero-size tridiagonal construction and solve.
- `NumericalIntegration::Instant` callable validation, finite input validation, interval handling, and normal integration.
- Legendre polynomial recurrence against analytic `P0` through `P3`.
- Determinant row-swap parity and input preservation.
- Owning `QinJiuShao::toFunction()` callable lifetime.
- SOR double relaxation factors, invalid omega, residual checks, and non-convergence.
- Multi-translation-unit link coverage remains in CTest.

M0.1 repaired:

- Upper/lower triangular copy construction and assignment now preserve base dimensions and storage.
- `TridiagonalMatrix(0)` no longer routes through `n - 1` allocation.
- `NumericalIntegration::Instant` checks the callable instead of assigning it to `nullptr`.
- Legendre recurrence uses floating-point coefficients.
- `determinant` applies row-swap parity.
- `QinJiuShao::toFunction()` returns an owning callable that captures coefficient state by value.
- `SOR(..., double omega, ...)` implements the standard update formula and validates `0 < omega < 2`; the legacy `int` overload remains.
- `LU_Decompose` error text now records the no-pivot limitation without mislabeling every zero-pivot failure as singular.

## M0.1 Actual Build/Test Result

Observed platform:

- OS shell: Windows PowerShell in `E:\Ducuments\codeForVScode\MatCal`.
- Compiler used by CMake: GNU C++ 15.1.0, `E:/minGW/mingw64/bin/c++.exe`.
- Build generator: `MinGW Makefiles`.

Commands run:

```powershell
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Debug
cmake --build build --config Debug
ctest --test-dir build --output-on-failure
cmake -S . -B build-release -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --config Release
ctest --test-dir build-release --output-on-failure
```

CTest result: 6/6 tests passed in both Debug and Release.

Sanitizers were not run in this M0.1 pass.

## Warning Baseline

See `docs/warning-baseline.md`. M0.1 intentionally fixes only warnings tied to the repaired bugs and safety contracts. Remaining warning families are preserved as documented legacy debt for M1 instead of being hidden by broad formatting or API churn.

## M6 Migration Audit

M6 reviewed `src/Iteration.hpp`, `src/Basics.hpp`, and `src/Insert.hpp` before changing migrated algorithms.

Reused legacy formulas and loops:

- `NewtonForEquations`: Newton system equation `J * delta = -F` and the finite-difference Jacobian idea.
- `Least_Square`: weighted normal-equation assembly, selected polynomial terms, and polynomial output shape.
- `LagrangeInsert`: Lagrange basis construction.
- `NewtonInsert_Quotient`: lower-triangular divided-difference recurrence.
- `NewtonInsert_Finite`: forward-difference table and finite-difference Newton basis.
- `Hermite`: squared Lagrange-basis Hermite value/derivative formula.

M6 repairs:

- `NewtonForEquations` no longer prints callable validation failures to stdout and no longer returns a last iterate as a pseudo successful root.
- `Least_Square` no longer hides rank-deficient or failed normal-equation solves behind placeholder success-like output.
- Remaining classic interpolation families now reject duplicate or non-finite nodes before denominator division or polynomial construction.
- `Insert.hpp` scalar `getDegree()` accessors no longer return `const int`.

New targets and tests:

- `MatCal::Nonlinear`
- `MatCal::LeastSquares`
- Extended `MatCal::Interpolation`
- Nonlinear, least-squares, remaining interpolation oracle tests.
- Legacy/new differential tests.
- Public-header self-contained tests.
- Multi-translation-unit tests.
- Nonlinear-only, LeastSquares-only, Interpolation-only, Legacy, PT-style, package, and add-subdirectory consumers.

## M7 Matrix Migration Audit

M7 reviewed the Matrix algorithms in `src/Matrix.hpp` before changing migrated paths.

Reused legacy formulas and loops:

- `Jacobi`: old-iterate stationary update over every component.
- `Gauss_Seidel`: lower entries from the current sweep, upper entries from the previous sweep.
- `SOR`: corrected relaxed Gauss-Seidel update already established in M0.1.
- `PowerMethod`: repeated matrix-vector products with vector normalization.
- `PowerMethod_reverse`: inverse iteration idea near a caller-provided shift, now using the dense partial-pivot solver instead of no-pivot legacy LU.
- `Derivative::pF_px` and `Derivative::dF_dx`: forward finite-difference partial derivative formula.

M7 repairs:

- Legacy Matrix-to-Linalg conversion now has one deep-copy adapter layer and does not rely on `dynamic_cast` to guess concrete matrix types.
- Legacy `solve_columnElimination` no longer owns a separate Gaussian-elimination implementation path; it delegates to `solve_dense_partial_pivot`.
- Stationary iterative solvers reject non-finite input, invalid omega, near-zero diagonals, dimension mismatches, and non-convergence with structured `SolverResult` state.
- Power and inverse-power solvers reject invalid initial vectors, non-finite values, singular shifted systems, and non-convergence without reporting pseudo success.
- Multivariable derivative helpers no longer mutate the caller-provided point while estimating partial derivatives.

New coverage:

- Legacy/Linalg adapter tests.
- Dense Jacobi, Gauss-Seidel, and SOR convergence, scale, and failure tests.
- Dense power and inverse-power eigen residual tests.
- Legacy/new differential tests for Matrix solvers, eigen solvers, and multivariable derivative helpers.
- Public-header self-contained and multi-translation-unit coverage for the new Linalg headers.

## M7.1 Numerical Contract Audit

M7.1 adds edge-case tests before making targeted fixes.

Confirmed bug:

- A finite stationary iterative system could produce an infinite initial residual when evaluating `Ax-b` for a finite nonzero initial guess. The solver did not immediately report `breakdown`; it continued into the iteration loop. The minimal regression case uses finite `DBL_MAX` matrix entries and initial guess `{1, 1}` so the first row residual sum overflows.

M7.1 repairs:

- Stationary solvers now return `SolverStatus::breakdown` immediately when initial residual metrics are non-finite, and they return no partial solution.
- Legacy `solve_columnElimination` now prevalidates every RHS entry before solving multi-RHS systems.
- The targeted legacy warning list was repaired without deleting public APIs: `QinJiuShao` signed/unsigned loop, deprecated private `sortByDegree()` call from `remove_fake()`, and focused-build `solve_Linear_System` unused warnings.

New coverage:

- Stationary solver zero-size, zero-matrix, 1x1, zero-RHS/nonzero-initial, invalid initial guess, omega-near-boundary, residual-overflow, input-preservation, and fixed-seed dense/direct differential tests.
- Eigen solver tiny/huge scale, near-zero initial vector, repeated dominant magnitude, near/exact shift, non-finite input, input-preservation, and fixed-seed symmetric residual tests.
- Legacy adapter empty roundtrip, non-finite rejection, multi-RHS direct solve, direct failure input preservation, determinant parity, and no-pivot LU failure tests.

## M8 Dense Pivoted LU Audit

M8 reviewed `solve_dense_partial_pivot`, legacy `solve_columnElimination`, legacy `determinant`, and legacy `LU_Decompose` before changing direct-solve paths.

Reused implementation ideas:

- Existing partial-pivot row search and row-swap elimination loop from the M1/M7 dense solver.
- Existing M1.1 scale-aware pivot tolerance, finite checks, residual acceptance contract, and `SolverMetrics` fields.
- Existing Legacy/Linalg deep-copy adapters for compatibility facades.

M8 repairs and refactors:

- `solve_dense_partial_pivot` no longer owns an independent full solve loop. It delegates to `factorize_dense_partial_pivot` and then solves one RHS.
- `PivotedLuFactorization` owns compact `LU`, row permutation, permutation sign, original matrix copy, and factorization metrics.
- Dense multi-RHS solve now factors once and solves all columns atomically.
- Legacy `solve_columnElimination` now uses one factorization for all RHS columns instead of re-running Gaussian elimination per column.
- Legacy `determinant` now delegates to pivoted LU and computes `permutation_sign * product(diag(U))` without modifying the input.
- Legacy `LU_Decompose` remains no-pivot LU because its old `LUresult` cannot represent the permutation in `P A = L U`.

New coverage:

- `PA = LU` reconstruction.
- Odd/even row-swap determinant signs.
- Singular, near-singular, non-finite, 0x0, 1x1, and scale `1e-20/1/1e20` cases.
- Multi-RHS reuse and operation-count checks.
- Input preservation, fixed-seed random solve cross-checks, public-header self-contained coverage, and multi-TU coverage.
