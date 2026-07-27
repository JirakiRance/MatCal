# MatCal::Linalg v0.1

M1 introduces `MatCal::Linalg` as a new 0.x development API. It is independent from the legacy `MatCal::Utils::Matrix` and `AbstractMatrix` hierarchy.

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

## Dense Reference Solver

`solve_dense_partial_pivot(const DenseMatrix&, const Vector&, const SolverOptions&)` is a small reference solver using partial-pivot Gaussian elimination.

It:

- Accepts only new `DenseMatrix` and `Vector`.
- Does not modify inputs.
- Returns `SolverResult` instead of printing.
- Checks dimensions, finite inputs, and solver options.
- Uses `absolute_tolerance + relative_tolerance * scale`, multiplied by `pivot_factor`, for pivot checks.
- Uses matrix infinity norm, clamped to at least `1`, as the dense scale.
- Computes infinity norm of `Ax - b`.
- Supports row-swap cases.
- Reports singular and near-singular pivots with structured status.

It is not a large-scale FEM solver and does not replace future Skyline, LDLT, CSR, or iterative solvers.
