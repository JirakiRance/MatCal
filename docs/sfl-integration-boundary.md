# SFL Integration Boundary

MatCal remains an independent Git repository. SFL may depend on MatCal, but MatCal must not depend on SFL.

Allowed dependency direction:

```text
SFL -> MatCal
```

## Integration Model

- SFL should pin MatCal through a Git submodule or another exact-version mechanism.
- MatCal should provide formal CMake exported targets, eventually including `MatCal::Linalg`.
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

## SFL Keeps

- FEM elements.
- Materials.
- Loads and boundary conditions.
- Mesh/Region ownership.
- Result recovery and interpretation.
- SFL diagnostics and source mapping.
