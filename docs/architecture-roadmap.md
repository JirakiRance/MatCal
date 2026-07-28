# Architecture Roadmap

M0 keeps the legacy surface intact and makes future change measurable. M0.1 adds targeted safety repairs and regression tests without a full rewrite. M1 introduces the first independent `MatCal::Linalg` 0.x API. M1.1 hardens its scale, finite-value, and result contracts. M2 adds general symmetric skyline storage and SPD LDLT. M2.1 optimizes and freezes that existing skyline baseline for later integration work. M3 starts the first Legacy modernization migration with a shared polynomial core and installable CMake package metadata. M4 migrates scalar roots plus linear and natural cubic spline interpolation. M5 migrates scalar calculus and ODE/RK4 cores. M6 migrates multivariable Newton systems, polynomial least squares, and remaining classic polynomial interpolation. M7 migrates legacy Matrix direct, stationary iterative, eigen, and multivariable derivative helper paths.

## Current Shape

MatCal now has compatibility and focused numerical targets:

- `MatCal::Legacy`: the header-only legacy library with broad public APIs in five headers.
- `MatCal::Linalg`: the new 0.x value-type linalg target under `include/MatCal/Linalg`.
- `MatCal::Polynomial`: the M3 value-type polynomial core extracted from `QinJiuShao`.
- `MatCal::Roots`: the M4 scalar root-finding core extracted from `Iteration.hpp`.
- `MatCal::Interpolation`: the M4 interpolation core extracted from `Insert.hpp`.
- `MatCal::Calculus`: the M5 finite-difference and quadrature core extracted from `Basics.hpp`.
- `MatCal::ODE`: the M5 Euler/RK4 core extracted from `Basics.hpp`.
- `MatCal::Nonlinear`: the M6 multivariable nonlinear-system core extracted from `Iteration.hpp`.
- `MatCal::LeastSquares`: the M6 polynomial least-squares core extracted from `Basics.hpp`.

`MatCal::Legacy` now depends on `MatCal::Polynomial`, `MatCal::Roots`, `MatCal::Interpolation`, `MatCal::Calculus`, `MatCal::ODE`, `MatCal::Nonlinear`, and `MatCal::LeastSquares` for migrated algorithms. `MatCal::Linalg` remains independent from Legacy.

## Direction

Target layering:

```text
legacy MatCal API
    -> compatibility facade
    -> internal linalg/numerics implementation
```

Public targets:

- `MatCal::Legacy`: current API surface.
- `MatCal::Linalg`: 0.x development API for independent vector/matrix storage and reference solvers.
- `MatCal::Polynomial`: 0.x development API for polynomial value operations.
- `MatCal::Roots`: 0.x development API for scalar nonlinear solves.
- `MatCal::Interpolation`: 0.x development API for interpolation.
- `MatCal::Calculus`: 0.x development API for scalar differentiation and integration.
- `MatCal::ODE`: 0.x development API for domain-neutral ODE stepping.
- `MatCal::Nonlinear`: 0.x development API for square nonlinear systems.
- `MatCal::LeastSquares`: 0.x development API for polynomial least-squares fitting.

Future public targets may include:

- `MatCal::Core`: scalar utilities and shared tolerance/status policies if they become useful outside linalg.

Do not describe `MatCal::Linalg` as ABI-stable or production-ready.

## Gradual Refactor Plan

1. Keep M0/M0.1 characterization and regression tests as the compatibility gate.
2. Maintain M1.1 `Vector`, `DenseMatrix`, `SolverOptions`, `SolverDiagnostic`, `SolverMetrics`, and `SolverResult` as the new baseline.
3. Maintain M2/M2.1 `SymmetricSkylineMatrix` and SPD `SkylineLdltFactorization` as general numerical facilities.
4. Keep the M3 legacy migration matrix current as each family moves from `legacy-only` to `migrated`.
5. Add pivoted LU factors as a new explicit API while preserving legacy no-pivot `LU_Decompose`.
6. Separate matrix storage from solver algorithms.
7. Replace stdout from core paths with structured status or caller-provided diagnostics.
8. Expand scale-aware tolerance helpers only when algorithms need them.
9. Forward legacy APIs to new internals one family at a time, behind compatibility tests.
10. Add CSR and advanced sparse iterative solver infrastructure after dense, skyline, and dense stationary behavior is stable.

## M1 Baseline

M1 builds on the M0.1 safety baseline rather than reopening the legacy headers wholesale:

- `Vector` and row-major `DenseMatrix` exist as independent value types.
- `SolverStatus`, `SolverOptions`, `SolverDiagnostic`, `SolverMetrics`, and `SolverResult` exist.
- `solve_dense_partial_pivot` is a reference dense solver with structured failure.
- Public linalg headers are tested for self-contained includes and multi-TU use.
- Legacy APIs are still separate and unchanged by M1.

M2 reuses:

- Scale-relative tolerance form.
- Explicit `matrix_scale`, `rhs_scale`, `solution_scale`, and residual metrics.
- `not_positive_definite` vs `breakdown` distinction for SPD LDLT.
- Machine-readable diagnostic `code`, `reason`, `phase`, `row`, and `column`.

M2.1 adds:

- Profile-limited skyline back substitution using factorization-owned column adjacency.
- Separate factorization and solve work counters in `SolverMetrics`.
- Performance contract tests based on deterministic storage/work counts, not wall-clock timing.
- The frozen 0.x tolerance composition `max(absolute_tolerance, relative_tolerance * scale)`.

M3 adds:

- `MatCal::Polynomial::Polynomial`.
- `MatCal::Legacy` delegation from `QinJiuShao` into the shared polynomial core.
- Installable CMake package files for `find_package(MatCal CONFIG)`.
- A migration matrix in `docs/legacy-migration-matrix.md`.

M4 adds:

- `MatCal::Roots` result/status/options contracts.
- Legacy delegation for `Bisection`, `Picard`, `Picard-Aitken`, `Newton`, downhill Newton, and secant variants.
- `MatCal::Interpolation::LinearInterpolator` and natural `CubicSpline`.
- Legacy delegation for `LinearInsert` and `CubicSpline`.
- PT-style source consumer and package/add-subdirectory consumer checks.

M5 adds:

- `MatCal::Calculus` result/status contracts for scalar finite differences and quadrature.
- Legacy delegation for `Derivative::dy_dx`, `Derivative::dy_dx_center`, `NumericalIntegration::Instant`, Newton-Cotes, composite Newton-Cotes, and Romberg.
- `MatCal::ODE` result/status contracts for simple Euler, improved Euler, and RK4 over vector states.
- Legacy delegation for `ODE::SimpleEuler`, `ODE::Euler`, `ODE::RungeKutta_44`, and PT-sensitive `Integrate::RK4::step/step2`.
- CMake export and consumer checks for `MatCal::Calculus` and `MatCal::ODE`.

M6 adds:

- `MatCal::Nonlinear` result/status/options contracts for square multivariable Newton solves.
- Legacy delegation for `NewtonForEquations::solve`.
- `MatCal::LeastSquares` result/status contracts for polynomial normal-equation fitting.
- Legacy delegation for all current `Least_Square::solve` overloads.
- `MatCal::Interpolation` polynomial interpolation helpers for Lagrange, divided-difference Newton, finite-difference Newton, and Hermite.
- Legacy delegation for `LagrangeInsert`, `NewtonInsert_Quotient`, `NewtonInsert_Finite`, and `Hermite`.
- CMake export and consumer checks for `MatCal::Nonlinear` and `MatCal::LeastSquares`.

M7 adds:

- Legacy/Linalg deep-copy conversion adapters owned on the Legacy side.
- `MatCal::Linalg::solve_jacobi`, `solve_gauss_seidel`, and `solve_sor`.
- `MatCal::Linalg::dominant_eigenpair` and `inverse_power_eigenpair`.
- Legacy delegation for `solve_columnElimination`, `Jacobi`, `Gauss_Seidel`, `SOR`, `PowerMethod`, and `PowerMethod_reverse`.
- `MatCal::Calculus::partial_difference` and `gradient`.
- Legacy delegation for `Derivative::pF_px` and `Derivative::dF_dx`.

## Future Capabilities

The architecture should eventually support:

- `Vector`
- `DenseMatrix`
- `SymmetricSkylineMatrix`
- `LDLT`
- `CSR`
- Advanced sparse/preconditioned iterative solvers

These are roadmap items, not M0 implementations.

## Non-Goals

- No FEM implementation.
- No SFL AST, diagnostics, regions, materials, elements, loads, or result IR.
- No mass rewrite of legacy headers.
- No ABI stability claim.
- No CSR, Skyline, LDLT, or FEM implementation in M1.
- No CSR, Skyline, LDLT, or FEM implementation in M1.1.
- No general indefinite symmetric factorization, Bunch-Kaufman, CSR, iterative method, FEM, or SFL integration in M2/M2.1.
