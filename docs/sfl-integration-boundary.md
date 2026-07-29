# SFL Integration Boundary

MatCal remains an independent Git repository. SFL may depend on MatCal, but MatCal must not depend on SFL.

Allowed dependency direction:

```text
SFL -> MatCal
```

## Integration Model

- SFL should pin MatCal through a Git submodule or another exact-version mechanism.
- MatCal now provides `MatCal::Linalg`, `MatCal::Polynomial`, `MatCal::Roots`, `MatCal::Interpolation`, `MatCal::Calculus`, `MatCal::ODE`, `MatCal::Nonlinear`, and `MatCal::LeastSquares` CMake targets as 0.x development APIs. SFL integration should wait for an explicit integration stage and pin an exact MatCal commit or tag.
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
- General scalar differentiation, quadrature, and ODE stepping over value-owned vectors.
- General multivariable nonlinear solving for numeric vectors.
- General polynomial least-squares fitting.
- General dense stationary iterative linear solves and dense power eigen solves.
- General dense partial-pivot LU factorization and multi-RHS dense solves.

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

M5 keeps calculus and ODE boundaries equally generic. `MatCal::ODE` accepts only numeric time, numeric state vectors, and numeric RHS callables. It does not accept PT mechanism types, SFL AST nodes, FEM element connectivity, materials, loads, constraints, stress recovery requests, CAE IR, diagnostics, or Result IR.

M6 keeps nonlinear and least-squares boundaries generic. `MatCal::Nonlinear` accepts numeric vectors, residual callables, optional Jacobian callables, and numeric options. `MatCal::LeastSquares` accepts numeric samples, weights, and polynomial term selections. Neither target accepts SFL AST, FEM element connectivity, constraints, mesh ownership, materials, loads, source diagnostics, CAE IR, or Result IR.

M7 keeps Matrix modernization generic. Legacy adapters convert only between legacy Matrix storage and `MatCal::Linalg` value types. Stationary iterative solvers and power eigen solvers accept only dense numeric matrices, vectors, numeric options, and numeric initial guesses. They do not accept SFL degrees of freedom, element connectivity, mesh constraints, materials, loads, source diagnostics, CAE IR, or Result IR.

M8 keeps direct dense factorization generic. `PivotedLuFactorization` accepts only numeric `DenseMatrix` input and returns generic solver metrics, diagnostics, row permutations, and solution matrices. It does not accept SFL degrees of freedom, assembled-element provenance, constraints, materials, loads, source diagnostics, CAE IR, or Result IR.

## SFL Keeps

- FEM elements.
- Materials.
- Loads and boundary conditions.
- Mesh/Region ownership.
- Result recovery and interpretation.
- SFL diagnostics and source mapping.
