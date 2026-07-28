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
| Interpolation | `LinearInsert`, `LagrangeInsert`, `NewtonInsert_*`, `CubicSpline`, `Hermite` | `src/Insert.hpp` | LinearInsert and natural CubicSpline migrated in M4; other interpolation families remain legacy-only | linear segment formula, binary interval search, natural spline tridiagonal equation, spline evaluation formula | strict finite/input checks, explicit extrapolation policy, no legacy matrix dependency for spline solve | `MatCal::Interpolation` | Yes for `LinearInsert` and `CubicSpline` | oracle, legacy differential, PT-style consumer | partially migrated |
| Scalar iteration | `Bisection`, `Picard`, `Newton`, `Secant` | `src/Iteration.hpp` | M4 migrated scalar loops and added structured failure states; bisection no longer accepts tiny scaled residual alone as endpoint success | bisection/Picard/Aitken/Newton/downhill/Secant formulas | finite checks, derivative/denominator zero, bracket checks, max-iteration status | `MatCal::Roots` | Yes for scalar root classes | oracle, failure, differential, header, multi-TU | migrated |
| Newton equations | `NewtonForEquations` | `src/Iteration.hpp` | Uses legacy Matrix, dynamic cast, fallback policy | finite-difference Jacobian idea | structured diagnostics, linalg direct solve, no stdout | none yet | No | legacy inventory only | planned |
| Differentiation | `Derivative::*` | `src/Basics.hpp` | Finite-difference helpers exist | central-difference formulas | finite checks, step validation | none yet | No | legacy inventory only | planned |
| Least squares | `Least_Square::solve` | `src/Basics.hpp` | Uses normal equations / legacy matrix | polynomial basis setup | conditioning/status, linalg solve | none yet | No | legacy inventory only | planned |
| Integration | `NumericalIntegration::*` | `src/Basics.hpp` | `Instant` fixed in M0.1; Newton-Cotes/Romberg still old | quadrature formulas | callable/finite contracts across all methods | none yet | No | Instant regression | characterized |
| ODE / RK4 | `ODE::SimpleEuler`, `Euler`, `RungeKutta_44` | `src/Basics.hpp` | RK4 stage formula exists; not yet extracted | stage formula and output shape | vector state API, dt/finite checks | none yet | No | legacy inventory only | planned |

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
