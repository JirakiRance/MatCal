# MatCal::Linalg v0.1

M1 introduces `MatCal::Linalg` as a new 0.x development API. M1.1 hardens its scale, finite-value, and `SolverResult` contracts. It is independent from the legacy `MatCal::Utils::Matrix` and `AbstractMatrix` hierarchy.

## Public Headers

- `include/MatCal/Linalg/Vector.hpp`
- `include/MatCal/Linalg/DenseMatrix.hpp`
- `include/MatCal/Linalg/SolverTypes.hpp`
- `include/MatCal/Linalg/DenseSolver.hpp`

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

## Dense Reference Solver

`solve_dense_partial_pivot(const DenseMatrix&, const Vector&, const SolverOptions&)` is a small reference solver using partial-pivot Gaussian elimination.

It:

- Accepts only new `DenseMatrix` and `Vector`.
- Does not modify inputs.
- Returns `SolverResult` instead of printing.
- Checks dimensions, finite inputs, and solver options.
- Uses `pivot_factor * (absolute_tolerance + relative_tolerance * matrix_scale)` for pivot checks.
- Uses maximum absolute matrix coefficient as `matrix_scale`; this value is not clamped to `1`.
- Defaults to `absolute_tolerance = 0`, `relative_tolerance = 1e-12`, `pivot_factor = 1`.
- Computes absolute residual infinity norm and relative residual.
- Supports row-swap cases.
- Reports singular and near-singular pivots with structured status.
- Returns `breakdown` for non-finite intermediate values during factorization, back substitution, or residual evaluation.
- Returns no partial-success solution on failure.

Scale invariance contract: multiplying a nonsingular finite system by ordinary finite factors such as `1e-20`, `1`, and `1e20` must not by itself change success into singular. This is covered by regression tests.

It is not a large-scale FEM solver and does not replace future Skyline, LDLT, CSR, or iterative solvers.
