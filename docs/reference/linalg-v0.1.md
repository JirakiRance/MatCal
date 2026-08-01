# MatCal::Linalg v0.1

M1 introduces `MatCal::Linalg` as a new 0.x development API. M1.1 hardens its scale, finite-value, and `SolverResult` contracts. M2 adds general symmetric skyline storage and SPD LDLT. M7 adds stationary iterative solvers and power eigen solvers. M8 adds reusable dense partial-pivot LU factorization. It is independent from the legacy `MatCal::Utils::Matrix` and `AbstractMatrix` hierarchy.

## Public Headers

- `include/MatCal/Linalg/Vector.hpp`
- `include/MatCal/Linalg/DenseMatrix.hpp`
- `include/MatCal/Linalg/SolverTypes.hpp`
- `include/MatCal/Linalg/DenseSolver.hpp`
- `include/MatCal/Linalg/IterativeSolvers.hpp`
- `include/MatCal/Linalg/EigenSolvers.hpp`
- `include/MatCal/Linalg/SymmetricSkylineMatrix.hpp`
- `include/MatCal/Linalg/SkylineLdlt.hpp`

Each public header is tested as self-contained and can be included without relying on include order.

## CMake Target

- Library target: `matcal_linalg`
- Alias target: `MatCal::Linalg`

`MatCal::Linalg` does not link to `MatCal::Legacy`, and `MatCal::Legacy` does not link to `MatCal::Linalg`. Future compatibility facades may bridge the two intentionally, but M1 keeps them parallel.

## Vector

`MatCal::Linalg::Vector` is an owning RAII value type backed by contiguous `std::vector<double>` storage.

Main contract:

- `std::size_t` sizes and indices.
- Default construction creates an empty vector.
- Size, fill-value, initializer-list, and explicit `std::vector<double>` construction.
- `operator[]` is unchecked and follows the standard container convention.
- `at()` checks bounds and throws `std::out_of_range`.
- `span()` returns borrowed views; callers must not outlive or reallocate the vector.
- `all_finite()` detects NaN/Inf.
- `dot`, `norm1`, stable scaled `norm2`, `normInf`, `axpy`, and `scale`.
- Size mismatches throw `std::invalid_argument`.

Ordinary vector operations do not return `SolverResult`. If finite inputs overflow during `dot`, `norm1`, `axpy`, or `scale`, the resulting double values follow IEEE behavior and may become Inf. `normInf()` reports Inf if present and NaN if any entry is NaN. Solvers must not continue with non-finite intermediates.

## DenseMatrix

`MatCal::Linalg::DenseMatrix` is an owning row-major dense matrix backed by one contiguous `std::vector<double>`. It does not inherit from `AbstractMatrix`.

Main contract:

- `std::size_t` rows and columns.
- Default construction creates a `0 x 0` matrix.
- `0 x N` and `N x 0` matrices are valid and store zero elements.
- Initializer-list construction rejects ragged rows.
- `operator()(row, col)` is unchecked.
- `at(row, col)` checks bounds and throws `std::out_of_range`.
- `row(i)` returns a borrowed row span; `N x 0` rows return empty spans.
- `checked_element_count()` rejects `rows * cols` overflow before allocation.
- `identity`, `fill`, `transpose`, matrix-vector multiply, matrix-matrix multiply, `normInf`, and `all_finite`.
- Dimension mismatches throw `std::invalid_argument`.

This dense container is a correctness baseline and small-problem utility. It is not the future storage for large FEM sparse systems.

`DenseMatrix::normInf()` is the conventional maximum absolute row sum and can return Inf if a finite row sum overflows. The dense solver does not use this as its pivot scale; it uses a safe maximum absolute coefficient scale so scale computation itself does not silently turn finite input into Inf.

## Dense Pivoted LU and Reference Solver

M8 adds reusable partial-pivot LU:

- `PivotedLuFactorization`
- `PivotedLuFactorizationResult`
- `MatrixSolverResult`
- `factorize_dense_partial_pivot`

The factorization owns the original matrix copy, compact `LU` data, row permutation, permutation sign, and factorization metrics. The contract is:

```text
P A = L U
```

where `L` is unit lower triangular in the strict lower part of `lu()`, `U` is the upper triangle of `lu()`, and `row_permutation()[i]` is the original row used for permuted row `i`.

`PivotedLuFactorization::solve(Vector)` and `solve(DenseMatrix)` reuse the same factorization. Matrix RHS solves are atomic: dimensions and finite RHS values are validated before computation, and numerical failure returns no partial solution matrix.

`solve_dense_partial_pivot(const DenseMatrix&, const Vector&, const SolverOptions&)` remains the one-shot reference solver, but it is now implemented as:

```text
factorize_dense_partial_pivot(A, options) -> factorization.solve(b, options)
```

It:

- Accepts only new `DenseMatrix` and `Vector`.
- Does not modify inputs.
- Returns `SolverResult` instead of printing.
- Checks dimensions, finite inputs, and solver options.
- Uses `pivot_factor * max(absolute_tolerance, relative_tolerance * matrix_scale)` for pivot checks.
- Uses maximum absolute matrix coefficient as `matrix_scale`; this value is not clamped to `1`.
- Defaults to `absolute_tolerance = 0`, `relative_tolerance = 1e-12`, `pivot_factor = 1`.
- Computes absolute residual infinity norm and relative residual.
- Supports row-swap cases.
- Reports singular and near-singular pivots with structured status.
- Returns `breakdown` for non-finite intermediate values during factorization, back substitution, or residual evaluation.
- Returns no partial-success solution on failure.
- Computes determinant through `permutation_sign * product(diag(U))`.

Scale invariance contract: multiplying a nonsingular finite system by ordinary finite factors such as `1e-20`, `1`, and `1e20` must not by itself change success into singular. This is covered by regression tests.

It is not a large-scale FEM solver and does not replace M2 skyline SPD LDLT, future CSR, or specialized iterative solvers.

## Stationary Iterative Solvers

M7 adds:

- `solve_jacobi`
- `solve_gauss_seidel`
- `solve_sor`

These functions reuse the legacy Jacobi, Gauss-Seidel, and corrected SOR update formulas over `DenseMatrix` and `Vector`. They return the existing `SolverResult` type and record iterations, residual metrics, matrix/RHS/solution scale, diagonal tolerance, and operation count.

Contracts:

- matrix must be square;
- RHS and optional initial guess must match size;
- inputs must be finite;
- diagonal entries must be above `pivot_factor * max(abs_tol, rel_tol * matrix_scale)`;
- `solve_sor` requires finite `0 < omega < 2`;
- success is based on the same absolute/relative residual acceptance style used by the dense solver;
- non-finite initial or per-iteration residual evaluation is reported as `SolverStatus::breakdown`;
- numerical failure and non-convergence do not return partial-success solutions.

## Power Eigen Solvers

M7 adds:

- `EigenStatus`
- `EigenOptions`
- `EigenDiagnostic`
- `EigenMetrics`
- `EigenResult`
- `dominant_eigenpair`
- `inverse_power_eigenpair`

The dominant power method reuses the legacy idea of repeatedly applying `A` and normalizing by the largest-magnitude component. The returned eigenvalue is computed by Rayleigh quotient and checked with residual `||Av - lambda v||_inf`.

The inverse power method applies shifted solves with `solve_dense_partial_pivot`; it does not keep the old no-pivot LU path. Singular shifted systems return `EigenStatus::singular_shift`.

These are dense 0.x reference algorithms, not a full eigensolver package. Repeated or clustered eigenvalues can converge slowly and may return `not_converged`.

## Symmetric Skyline and SPD LDLT

M2 adds `SymmetricSkylineMatrix` and `SkylineLdltFactorization`; M2.1 freezes the first integration-ready SPD skyline contract.

`SymmetricSkylineMatrix` is a general numerical symmetric matrix container with row-oriented lower skyline storage. It supports explicit first-column profiles and symmetric nonzero position pairs. It does not accept FEM element connectivity or any SFL type.

`factorize_skyline_ldlt()` computes an unpivoted SPD-only LDLT factorization. `solve_skyline_ldlt()` is a one-shot helper; callers that need multiple right-hand sides should call `factorize_skyline_ldlt()` once and reuse the returned `SkylineLdltFactorization`.

M2.1 keeps the same public classes and functions. It optimizes existing back substitution to use column profile adjacency rather than scanning all later rows, and it extends `SolverMetrics` with `factorization_operation_count` and `solve_operation_count` while keeping `operation_count` as a total/reference counter.

See `docs/reference/skyline-ldlt.md` for the full storage and numerical contract.

## M8.1 API Audit Notes

The `MatCal::Linalg` public surface for `v0.3.0-alpha.1` is still a 0.x source API:

- `PivotedLuFactorization` is copyable and movable by value because it owns standard-library storage only.
- `PivotedLuFactorization::lu()`, `row_permutation()`, and `permutation_sign()` expose the numerical factorization contract without borrowing from input matrices.
- `PivotedLuFactorization::solve(...)` does not mutate the factorization, so multiple RHS solves reuse the same data.
- Failure results do not contain partial solutions.
- Dense, skyline, iterative, and eigen APIs use structured diagnostics and do not depend on SFL diagnostics.
