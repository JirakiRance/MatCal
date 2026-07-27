# SFL Integration Boundary

MatCal remains an independent Git repository. SFL may depend on MatCal, but MatCal must not depend on SFL.

Allowed dependency direction:

```text
SFL -> MatCal
```

## Integration Model

- SFL should pin MatCal through a Git submodule or another exact-version mechanism.
- MatCal now provides a `MatCal::Linalg` CMake target as a 0.x development API. SFL integration should wait for an explicit integration stage and pin an exact MatCal commit or tag.
- SFL should consume MatCal with `target_link_libraries`.
- SFL should not copy MatCal source files into its own repository.
- SFL should not maintain a forked MatCal subtree as its internal implementation.

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

## M1 Boundary

M1 is not an SFL integration release. It adds `Vector`, `DenseMatrix`, structured solver types, and a dense reference solver only.

Skyline and LDLT are deferred to M2 or later. When added, they must remain general matrix/solver facilities, not SFL-owned mesh, element, material, load, or result-recovery logic.

## SFL Keeps

- FEM elements.
- Materials.
- Loads and boundary conditions.
- Mesh/Region ownership.
- Result recovery and interpretation.
- SFL diagnostics and source mapping.
