# MatCal v0.3.0-alpha.2 API Manifest

This manifest records the public source API expected for the `v0.3.0-alpha.2` release candidate. It is intended to be source-compatible with `v0.3.0-alpha.1`.

It is a compatibility audit baseline, not an ABI stability promise.

## Version and ABI

- CMake project version: `0.3.0`
- Proposed release tag: `v0.3.0-alpha.2`
- Previous tag: `v0.3.0-alpha.1`
- ABI stability: not promised
- Default integration model: pinned source build through CMake targets

## Public API Delta From alpha.1

No public numerical API change is intended in alpha.2.

The expected differences from alpha.1 are release engineering only:

- improved Linux CI diagnostics;
- alpha API manifest documentation on main;
- explicit release-asset notes.

## CMake Targets

- `MatCal::Legacy`
- `MatCal::Linalg`
- `MatCal::Polynomial`
- `MatCal::Roots`
- `MatCal::Interpolation`
- `MatCal::Calculus`
- `MatCal::ODE`
- `MatCal::Nonlinear`
- `MatCal::LeastSquares`

`MatCal::Legacy` exports dependencies on migrated header-only cores. `MatCal::Linalg` does not depend on `MatCal::Legacy`.

## Installed Public Headers

New public headers under `include/MatCal`:

- `MatCal/Calculus/Calculus.hpp`
- `MatCal/Interpolation/CubicSpline.hpp`
- `MatCal/Interpolation/LinearInterpolator.hpp`
- `MatCal/Interpolation/PolynomialInterpolation.hpp`
- `MatCal/LeastSquares/LeastSquares.hpp`
- `MatCal/Linalg/DenseMatrix.hpp`
- `MatCal/Linalg/DenseSolver.hpp`
- `MatCal/Linalg/EigenSolvers.hpp`
- `MatCal/Linalg/IterativeSolvers.hpp`
- `MatCal/Linalg/SkylineLdlt.hpp`
- `MatCal/Linalg/SolverTypes.hpp`
- `MatCal/Linalg/SymmetricSkylineMatrix.hpp`
- `MatCal/Linalg/Vector.hpp`
- `MatCal/Nonlinear/Nonlinear.hpp`
- `MatCal/ODE/ODE.hpp`
- `MatCal/Polynomial/Polynomial.hpp`
- `MatCal/Roots/Roots.hpp`

Legacy root headers installed for source compatibility:

- `Basics.hpp`
- `Insert.hpp`
- `Iteration.hpp`
- `Matrix.hpp`
- `QinJiuShao.hpp`

## MatCal::Linalg

Primary public types:

- `Vector`
- `DenseMatrix`
- `SolverStatus`
- `SolverOptions`
- `SolverDiagnostic`
- `SolverMetrics`
- `SolverResult`
- `MatrixSolverResult`
- `PivotedLuFactorization`
- `PivotedLuFactorizationResult`
- `SymmetricSkylineMatrix`
- `SkylineLdltFactorization`
- `SkylineLdltFactorizationResult`
- `EigenStatus`
- `EigenOptions`
- `EigenDiagnostic`
- `EigenMetrics`
- `EigenResult`

Primary public functions:

- `residual_norm_inf`
- `factorize_dense_partial_pivot`
- `solve_dense_partial_pivot`
- `factorize_skyline_ldlt`
- `solve_skyline_ldlt`
- `solve_jacobi`
- `solve_gauss_seidel`
- `solve_sor`
- `dominant_eigenpair`
- `inverse_power_eigenpair`
- `to_string(SolverStatus)`
- `to_string(EigenStatus)`

Default `SolverOptions`:

- `absolute_tolerance = 0`
- `relative_tolerance = 1e-12`
- `pivot_factor = 1.0`
- `max_iterations = 1000`

Solver status vocabulary:

- `success`
- `invalid_input`
- `dimension_mismatch`
- `non_finite_input`
- `singular`
- `not_positive_definite`
- `breakdown`
- `not_converged`

## Other Public Targets

The following targets retain the alpha.1 source API:

- `MatCal::Polynomial`: `Polynomial`
- `MatCal::Roots`: scalar root result/options/status types and bisection, Picard, Newton, downhill Newton, and secant functions
- `MatCal::Interpolation`: `LinearInterpolator`, `CubicSpline`, Lagrange/Newton/Hermite polynomial interpolation helpers
- `MatCal::Calculus`: scalar finite differences, gradient helpers, Newton-Cotes, composite Newton-Cotes, and Romberg
- `MatCal::ODE`: Euler, improved Euler, RK4 step and trajectory helpers
- `MatCal::Nonlinear`: Newton system solve with analytic or finite-difference Jacobian
- `MatCal::LeastSquares`: polynomial normal-equation fitting helpers

## MatCal::Legacy

Legacy source compatibility remains centered on:

- `MatCal::Utils`
- `MatCal::Algorithm::Matrix`
- `MatCal::Algorithm::Insert`
- `MatCal::Algorithm::Iteration`
- `MatCal::Algorithm::Basics`
- `MatCal::Algorithm::Basics::Integrate`

Public legacy headers, class names, constructor forms, and main function names are retained. Migrated facades forward to modern cores where compatibility allows. `LU_Decompose` remains legacy no-pivot LU because its result cannot represent the permutation in `P A = L U`.

## Release Asset Notes

- `.gitmodules`: absent; no submodules are required.
- License at the published tag snapshot: no `LICENSE` file was present.
- License on current `main` after D1 maintenance cleanup: MIT.
- The published tag source archive contains CMake/package files, `include/`,
  `src/`, `tests/`, `docs/`, `.github/workflows/ci.yml`, `README.md`, and
  `.gitignore`.
- Current `main` after D1 also contains `LICENSE` and installs it with the CMake
  package.
- Build directories and generated install trees must remain untracked.
