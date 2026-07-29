# MatCal v0.3.0-alpha.1 Release Notes

`v0.3.0-alpha.1` is the planned release candidate for the migrated MatCal 0.x source API after M8. It must be tagged only after the release-preparation commit is pushed and the remote Linux/Windows CI jobs are green.

## Version

- CMake project version: `0.3.0`
- Proposed Git tag: `v0.3.0-alpha.1`
- Previous integration baseline: `v0.2.0-alpha.1` for the M2.1 Skyline/LDLT baseline

MatCal remains a 0.x source API. ABI stability is not promised; downstream projects should rebuild from a pinned commit or tag.

## Public Targets

- `MatCal::Legacy`
- `MatCal::Linalg`
- `MatCal::Polynomial`
- `MatCal::Roots`
- `MatCal::Interpolation`
- `MatCal::Calculus`
- `MatCal::ODE`
- `MatCal::Nonlinear`
- `MatCal::LeastSquares`

## Highlights Since v0.2.0-alpha.1

- `MatCal::Polynomial` value core and QinJiuShao delegation.
- `MatCal::Roots` scalar root solving with structured results.
- `MatCal::Interpolation` linear, natural cubic spline, Lagrange, Newton, finite-difference Newton, and Hermite interpolation cores.
- `MatCal::Calculus` finite differences and quadrature.
- `MatCal::ODE` Euler and RK4 stepping/integration over value-owned numeric states.
- `MatCal::Nonlinear` square Newton systems with analytic or finite-difference Jacobian.
- `MatCal::LeastSquares` polynomial normal-equation fitting with structured status.
- `MatCal::Linalg` stationary iterative solvers, power/inverse-power eigen solvers, and reusable dense partial-pivot LU.
- Legacy facades now delegate many formerly duplicated algorithms to the new cores while preserving source-level legacy entry points.

## M8 Direct Solver Changes

- `PivotedLuFactorization` owns compact `LU`, row permutation, permutation sign, original matrix copy, and metrics.
- `factorize_dense_partial_pivot` exposes reusable partial-pivot factorization.
- `solve_dense_partial_pivot` is a one-shot facade over factorization plus solve.
- Multi-RHS dense solves factor once and are atomic.
- Legacy `solve_columnElimination` now reuses one factorization for all RHS columns.
- Legacy `determinant` uses `permutation_sign * product(diag(U))`.
- Legacy `LU_Decompose` remains no-pivot because its old result type cannot represent `P` in `P A = L U`.

## Build and Package

- Standard CMake install and `find_package(MatCal CONFIG REQUIRED)` are supported.
- Package consumers may request `find_package(MatCal 0.3 CONFIG REQUIRED)`.
- The default downstream integration path is a pinned source build through CMake targets, not copying MatCal sources and not using precompiled binaries by default.
- GitHub Actions workflow is added for Linux GCC and Windows MSVC Debug/Release builds, CTest, install, package consumer, add-subdirectory consumer, and SFL-compatible consumer checks.
- Remote CI must run before claiming the release candidate is verified.

## SFL Upgrade Note

SFL may evaluate upgrading from `v0.2.0-alpha.1` to `v0.3.0-alpha.1` by pinning the new tag or exact commit and linking MatCal targets. MatCal still accepts only generic numerical inputs and returns generic numerical results. It does not contain SFL AST, FEM element/material/load concepts, SFL diagnostics, CAE IR, or Result IR.

## PT Compatibility Note

Legacy headers, namespaces, constructors, and PT-sensitive APIs such as interpolation classes and RK4 `step`/`step2` remain available through `MatCal::Legacy`. Source consumers should recompile when upgrading; ABI stability for old compiled binaries is not promised.

## Known Limits

- No ABI stability promise.
- No CSR, QR, SVD, Bunch-Kaufman, or general indefinite symmetric factorization.
- Skyline LDLT remains SPD-only and unpivoted.
- Dense iterative and eigen solvers are reference algorithms, not large-scale production solvers.
- Least squares currently uses normal equations.
- `LU_Decompose` remains legacy no-pivot LU.
