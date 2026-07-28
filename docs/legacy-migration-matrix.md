# Legacy Migration Matrix

Stage M3 starts treating Legacy as a source of reviewed algorithms rather than as a sealed compatibility island. This matrix is intentionally operational: each row records where the old behavior lives, what can be reused, and whether a new core already owns the algorithm.

## Status Legend

- `characterized`: covered by legacy characterization or regression tests.
- `migrated`: new core exists and legacy delegates at least the shared algorithm.
- `planned`: reviewed enough to schedule, not yet migrated.
- `legacy-only`: still independent old implementation.

## Matrix

| Area | Legacy API | Algorithm Location | Correctness / Bugs | Reusable Parts | Needs Repair | New Core | Legacy Delegates | Tests | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Dense matrix | `MatCal::Utils::Matrix` | `src/Matrix.hpp` | Basic value storage characterized; some constness and dynamic-cast risks remain | Constructors, indexing behavior, arithmetic characterization | ownership, const-correctness, conversion checks | `MatCal::Linalg::DenseMatrix` | No | characterization, linalg oracle | characterized |
| Triangular matrix | `UpperTriangularMatrix`, `LowerTriangularMatrix` | `src/Matrix.hpp` | M0.1 fixed copy/base dimensions | storage and solve behavior | future facade to value storage | none yet | No | triangular regression | characterized |
| Tridiagonal matrix | `TridiagonalMatrix` | `src/Matrix.hpp` | M0.1 fixed zero-size construction | Thomas/chase behavior | structured status | none yet | No | zero-size regression | characterized |
| Gaussian elimination | `solve_columnElimination` | `src/Matrix.hpp` | Partial pivot path exists; dimension/pivot diagnostics incomplete | row-pivot elimination | structured result, finite checks | `MatCal::Linalg::solve_dense_partial_pivot` | No | legacy characterization, linalg oracle | characterized |
| LU | `LU_Decompose` | `src/Matrix.hpp` | No-pivot LU only; zero leading pivot can mean pivoting required | Doolittle loop for no-pivot cases | new pivoted API later | none yet | No | M0/M0.1 tests | characterized |
| Jacobi / Gauss-Seidel / SOR | `Jacobi`, `Gauss_Seidel`, `SOR` | `src/Matrix.hpp` | SOR omega/update fixed in M0.1 | iteration skeletons | common iterative result contract | none yet | No | SOR regression | characterized |
| Eigenvalue | power / inverse power methods | `src/Matrix.hpp` | Not yet migrated in M3 | legacy iteration idea | status, finite checks, shift/inverse diagnostics | none yet | No | legacy inventory only | planned |
| Polynomial | `QinJiuShao`, `QinJiuShaoNode` | `src/QinJiuShao.hpp` | Horner, arithmetic, derivative, integral are mathematically reusable; old node constructor printed before throwing | coefficient cleanup, output facade, constructor compatibility | finite checks, owning callable, central algorithm ownership | `MatCal::Polynomial::Polynomial` | Yes for evaluate, callable, add/subtract/multiply, scalar ops, derivative, integral, definite integral | new oracle, differential, negative, header, multi-TU | migrated |
| Interpolation | `LinearInsert`, `LagrangeInsert`, `NewtonInsert_*`, `CubicSpline`, `Hermite` | `src/Insert.hpp` | LinearInsert and natural CubicSpline migrated in M4; remaining classic polynomial interpolation migrated in M6 | linear segment formula, binary interval search, natural spline tridiagonal equation, Lagrange basis, divided differences, finite differences, Hermite basis | no barycentric/rational/multidimensional interpolation yet | `MatCal::Interpolation` | Yes | oracle, legacy differential, PT-style consumer, header, multi-TU | migrated |
| Scalar iteration | `Bisection`, `Picard`, `Newton`, `Secant` | `src/Iteration.hpp` | M4 migrated scalar loops and added structured failure states; bisection no longer accepts tiny scaled residual alone as endpoint success | bisection/Picard/Aitken/Newton/downhill/Secant formulas | finite checks, derivative/denominator zero, bracket checks, max-iteration status | `MatCal::Roots` | Yes for scalar root classes | oracle, failure, differential, header, multi-TU | migrated |
| Newton equations | `NewtonForEquations` | `src/Iteration.hpp` | M6 migrates Newton equation and finite-difference Jacobian idea; old pseudo-success failure paths removed | `J * delta = -F`, finite-difference Jacobian | line search/damping not yet implemented | `MatCal::Nonlinear` | Yes | oracle, failure states, legacy differential, header, multi-TU | migrated |
| Differentiation | `Derivative::*` | `src/Basics.hpp` | `dy_dx` and `dy_dx_center` migrated in M5; multivariable helpers remain legacy-only | forward and centered finite-difference formulas | finite checks, step validation | `MatCal::Calculus` | Yes for scalar `dy_dx` and `dy_dx_center` | oracle, invalid input, legacy differential | partially migrated |
| Least squares | `Least_Square::solve` | `src/Basics.hpp` | M6 migrates weighted and selected-term normal equations; rank deficiency no longer returns pseudo-success | Vandermonde/normal-equation construction, weights, selected terms, polynomial output | QR/SVD future API for conditioning | `MatCal::LeastSquares` | Yes | exact/noisy/weighted/selected/rank tests, legacy differential | migrated |
| Integration | `NumericalIntegration::*` | `src/Basics.hpp` | M5 migrates Instant, Newton-Cotes, Composite Newton-Cotes, and Romberg; discrete-data overload forwards through LinearInsert facade | left rectangle, closed Newton-Cotes weights, composite loop, Romberg trapezoid/Richardson loop | signed interval contract, non-finite checks, non-convergence status | `MatCal::Calculus` | Yes for callable overloads | oracle, invalid input, legacy differential | migrated |
| ODE / RK4 | `ODE::SimpleEuler`, `Euler`, `RungeKutta_44`, `Integrate::RK4::step/step2` | `src/Basics.hpp` | M5 migrates simple Euler, improved Euler, classic RK4, and PT-style RK4 helpers | Euler formulas, RK4 four-stage formula, legacy trajectory row shape | vector-state API, finite checks, RHS size checks, structured diagnostics | `MatCal::ODE` | Yes | oracle, PT-style differential, legacy differential | migrated |

## M3 Completed Migration

`QinJiuShao` now adapts to `MatCal::Polynomial::Polynomial` for the core polynomial algorithms:

- Horner evaluation;
- owning callable;
- add/subtract/multiply;
- scalar multiply/divide;
- derivative;
- indefinite and definite integral.

The old string formatting, `show()`, sparse descending node list, constructors, and public class names remain as the compatibility facade.

## M4 Completed Migration

`MatCal::Roots` now owns scalar root-finding loops for bisection, Picard, Picard-Aitken, Newton, downhill Newton, and secant variants. Legacy scalar classes delegate to the new core.

`MatCal::Interpolation` now owns linear interpolation and natural cubic spline interpolation. Legacy `LinearInsert` and `CubicSpline` delegate to the new core while preserving PT-sensitive old names, constructors, `calculate`, and getters.

## M5 Completed Migration

`MatCal::Calculus` now owns scalar forward/central finite differences plus left-rectangle, Newton-Cotes, composite Newton-Cotes, and Romberg integration. Legacy scalar derivative and callable integration APIs delegate to this core.

`MatCal::ODE` now owns simple Euler, improved Euler, and classic RK4 stepping/integration over value-owned state vectors. Legacy `ODE::*` table APIs and PT-sensitive `Integrate::RK4::step/step2` delegate to the shared RK4/Euler core.

## M6 Completed Migration

`MatCal::Nonlinear` now owns the multivariable Newton loop for square systems. Legacy `NewtonForEquations::solve` delegates to the finite-difference path while preserving its public result struct.

`MatCal::LeastSquares` now owns polynomial normal-equation fitting for complete and selected terms. Legacy `Least_Square::solve` overloads delegate to this core while converting results back to legacy `Matrix` and `QinJiuShao` objects.

`MatCal::Interpolation` now owns Lagrange, divided-difference Newton, finite-difference Newton, and Hermite polynomial interpolation. The corresponding legacy classes delegate to this core while keeping their constructors, `calculate`, `reconstruct`, `getPoly`, and table getters.
