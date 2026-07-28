# SFL Integration Boundary

MatCal remains an independent Git repository. SFL may depend on MatCal, but MatCal must not depend on SFL.

Allowed dependency direction:

```text
SFL -> MatCal
```

## Integration Model

- SFL should pin MatCal through a Git submodule or another exact-version mechanism.
- MatCal now provides `MatCal::Linalg`, `MatCal::Polynomial`, `MatCal::Roots`, and `MatCal::Interpolation` CMake targets as 0.x development APIs. SFL integration should wait for an explicit integration stage and pin an exact MatCal commit or tag.
- SFL should consume MatCal with `target_link_libraries` and the needed MatCal targets.
- SFL should not copy MatCal source files into its own repository.
- SFL should not maintain a forked MatCal subtree as its internal implementation.
- SFL should not use a precompiled MatCal binary as the default integration path; the default should be a pinned source version built through CMake.

## MatCal Must Not Contain

- SFL AST
- CAE IR
- SFL diagnostics
- Mesh or Region concepts
- Material, element, or load semantics
- FEM-specific business logic
- Result IR

## MatCal May Provide

- General matrices and vectors.
- General solver options.
- General solver results.
- Numerical algorithms with domain-neutral contracts.
- General scalar root solving and interpolation.

## M1 Boundary

M1/M1.1/M2/M2.1 is not an SFL integration release. It adds `Vector`, `DenseMatrix`, `SymmetricSkylineMatrix`, structured solver types, diagnostics, metrics, a dense reference solver, and an SPD skyline LDLT baseline only.

Skyline and LDLT remain general matrix/solver facilities, not SFL-owned mesh, element, material, load, or result-recovery logic.

M1.1 diagnostics remain MatCal-owned and domain-neutral. They are not SFL diagnostics and do not carry SFL AST, source mapping, mesh, region, material, element, load, or result IR fields.

M2 skyline profiles are built from first-column profiles or symmetric nonzero position pairs. SFL may later translate its own assembly structure into these generic numerical inputs, but MatCal does not accept SFL element connectivity or own the assembly semantics.

M2.1 freezes the generic capabilities needed for a later integration stage:

- create a skyline matrix from a first-column profile;
- accumulate entries with `add(row, column, value)`;
- multiply by `Vector`;
- factorize once and solve multiple right-hand sides;
- override `SolverOptions`;
- read matrix scale, pivot tolerance, minimum pivot, residuals, and work counters from `SolverMetrics`;
- receive structured non-success results without partial solutions for numerical failures.

## SFL Keeps

- FEM elements.
- Materials.
- Loads and boundary conditions.
- Mesh/Region ownership.
- Result recovery and interpretation.
- SFL diagnostics and source mapping.
