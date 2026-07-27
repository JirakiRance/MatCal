# Skyline Storage and SPD LDLT

M2 adds a general numerical skyline storage type and a baseline SPD LDLT solver under `MatCal::Linalg`.

## Public API

- `SymmetricSkylineMatrix`
- `SkylineLdltFactorization`
- `SkylineLdltFactorizationResult`
- `factorize_skyline_ldlt(...)`
- `solve_skyline_ldlt(...)`

These APIs do not depend on `MatCal::Legacy`, `AbstractMatrix`, SFL, FEM, mesh, node, element, degree-of-freedom, load, material, CAE IR, or Result IR types.

## Skyline Storage

`SymmetricSkylineMatrix` stores one real symmetric matrix using row-oriented lower skyline storage.

For row `i`, `first_columns()[i]` gives the first stored column in that row. Columns from `first_columns()[i]` through `i` are stored contiguously. Access is symmetric, so `get(r, c)` and `get(c, r)` observe the same stored value.

Construction supports:

- explicit first-column profile;
- symmetric nonzero position pairs.

Position pairs describe numerical matrix structure only. They are not FEM connectivity and do not encode element, node, region, or constraint semantics.

Contract:

- `get()` outside the profile returns structural zero.
- `set()` and `add()` outside the profile throw `std::out_of_range`.
- profile gaps inside the skyline are stored and default to zero.
- `multiply()` is matrix-vector multiply over the symmetric stored profile.
- `storage_size()` is `sum(i - first_columns[i] + 1)`.
- all index and storage-size computations are checked for invalid profiles and overflow.
- `matrix_scale()` is the maximum absolute stored coefficient.
- `all_finite()` validates stored values.

Storage complexity is `O(storage_size)`. Matrix-vector multiply is `O(storage_size)`.

## SPD LDLT

`factorize_skyline_ldlt()` computes an unpivoted LDLT factorization for real symmetric positive definite matrices:

```text
A = L D L^T
```

where `L` is unit lower triangular and stored in the same skyline profile, and `D` is a positive diagonal vector.

Contract:

- input matrix is not modified;
- factorization owns its factor data;
- solve is separated from factorization, so one factorization can solve multiple right-hand sides;
- the implementation does not convert to `DenseMatrix` for factorization or solve;
- only SPD matrices are supported in M2.

Non-positive or too-small pivots return `SolverStatus::not_positive_definite` with diagnostic code `non_positive_pivot`. Non-finite intermediate values return `SolverStatus::breakdown`. Non-finite inputs return `SolverStatus::non_finite_input`.

Factorization work is roughly:

```text
O(sum(profile_width_i^2))
```

Solve is currently a reference implementation. Forward substitution and factor matvec are skyline-profile operations; back substitution scans later rows and is not yet optimized for very large systems.

## Metrics and Residual

Skyline LDLT reuses `SolverOptions`, `SolverResult`, `SolverDiagnostic`, and `SolverMetrics`.

Metrics include:

- `matrix_scale`
- `rhs_scale`
- `solution_scale`
- `pivot_tolerance_used`
- `minimum_abs_pivot`
- `absolute_residual_norm`
- `relative_residual_norm`
- `residual_acceptance_tolerance`
- `operation_count`

Residual uses the M1.1 contract:

```text
abs_res = ||Ax-b||_inf
rel_res = abs_res / (matrix_scale * solution_scale + rhs_scale)
```

Success requires:

```text
abs_res <= absolute_tolerance + relative_tolerance * (matrix_scale * solution_scale + rhs_scale)
```

## Deferred Work

M2 does not implement:

- indefinite symmetric factorization;
- Bunch-Kaufman pivoting;
- CSR;
- iterative methods;
- FEM assembly;
- SFL integration.
