# MatCal v0.3.0-alpha.1 API Manifest

This manifest freezes the public source API snapshot for the `v0.3.0-alpha.1` release candidate. It is a compatibility audit baseline, not an ABI stability promise.

## Version and ABI

- CMake project version: `0.3.0`
- Proposed release tag: `v0.3.0-alpha.1`
- Candidate commit: `31f493789a5b2f191aef59c828ec7562a7f9412f`
- ABI stability: not promised
- Default integration model: pinned source build through CMake targets

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

## MatCal::Polynomial

Primary public type:

- `Polynomial`

The polynomial core owns dense coefficient storage and supports value evaluation, arithmetic, derivative, integral, definite integral, and owning callable conversion. Very high sparse degrees are guarded before dense allocation.

## MatCal::Roots

Primary public types:

- `RootStatus`
- `RootReason`
- `RootPhase`
- `RootOptions`
- `RootDiagnostic`
- `RootMetrics`
- `RootResult`

Primary public functions:

- `solve_bisection`
- `solve_picard`
- `solve_picard_aitken`
- `solve_newton`
- `solve_downhill_newton`
- `solve_secant_two_point`
- `solve_secant_one_point`

## MatCal::Interpolation

Primary public types:

- `ExtrapolationPolicy`
- `LinearInterpolator`
- `CubicSpline`
- `DividedDifferenceResult`

Primary public functions:

- `interpolate_lagrange`
- `interpolate_newton_divided`
- `interpolate_newton_finite`
- `interpolate_hermite`

## MatCal::Calculus

Primary public types:

- `CalculusStatus`
- `CalculusReason`
- `CalculusPhase`
- `CalculusDiagnostic`
- `DerivativeResult`
- `GradientResult`
- `IntegrationOptions`
- `IntegrationMetrics`
- `IntegrationResult`

Primary public functions:

- `forward_difference`
- `central_difference`
- `partial_difference`
- `gradient`
- `integrate_instant`
- `integrate_newton_cotes`
- `integrate_composite_newton_cotes`
- `integrate_romberg`

## MatCal::ODE

Primary public types:

- `OdeStatus`
- `OdeReason`
- `OdePhase`
- `OdeDiagnostic`
- `OdeOptions`
- `OdeMetrics`
- `OdeStepResult`
- `OdeTrajectoryResult`

Primary public functions:

- `rk4_step`
- `euler_step`
- `improved_euler_step`
- `integrate_rk4`
- `integrate_euler`

## MatCal::Nonlinear

Primary public types:

- `NonlinearStatus`
- `NonlinearReason`
- `NonlinearPhase`
- `NonlinearOptions`
- `NonlinearDiagnostic`
- `NonlinearMetrics`
- `NonlinearResult`

Primary public functions:

- `solve_newton_system`
- `finite_difference_jacobian`
- `solve_newton_system_finite_difference`

## MatCal::LeastSquares

Primary public types:

- `LeastSquaresStatus`
- `LeastSquaresReason`
- `LeastSquaresDiagnostic`
- `LeastSquaresMetrics`
- `LeastSquaresResult`

Primary public functions:

- `fit_polynomial`
- `fit_polynomial_degree`
- `fit_polynomial_selected`

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
- `LICENSE`, `NOTICE`, and `THIRD_PARTY` files: absent in this repository snapshot.
- Source archive should contain CMake/package files, `include/`, `src/`, `tests/`, `docs/`, `.github/workflows/ci.yml`, `README.md`, and `.gitignore`.
- Build directories and generated install trees must remain untracked.
