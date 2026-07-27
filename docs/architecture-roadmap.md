# Architecture Roadmap

M0 keeps the legacy surface intact and makes future change measurable. It intentionally does not perform a full rewrite.

## Current Shape

MatCal is a header-only legacy library with broad public APIs in five headers. Storage types, algorithms, printing, exceptions, and solver policy are tightly coupled.

## Direction

Target layering:

```text
legacy MatCal API
    -> compatibility facade
    -> internal linalg/numerics implementation
```

Future public targets may include:

- `MatCal::Legacy`: current API surface.
- `MatCal::Core`: scalar utilities, status/result types, tolerance policy.
- `MatCal::Linalg`: matrix/vector storage and linear solvers.

Do not publish `MatCal::Core` or `MatCal::Linalg` as stable until the contracts exist.

## Gradual Refactor Plan

1. Freeze M0 tests and docs.
2. Introduce `Vector`, `DenseMatrix`, and `SolverStatus/SolverResult` internally.
3. Separate matrix storage from solver algorithms.
4. Replace stdout from core paths with structured status or caller-provided diagnostics.
5. Add scale-aware tolerance helpers.
6. Forward legacy APIs to new internals one family at a time.
7. Add CSR and iterative solver infrastructure after dense behavior is stable.
8. Add `SymmetricSkylineMatrix` and `LDLT` for future structural mechanics needs, without adding FEM semantics.

## Future Capabilities

The architecture should eventually support:

- `Vector`
- `DenseMatrix`
- `SymmetricSkylineMatrix`
- `LDLT`
- `CSR`
- Iterative solvers

These are roadmap items, not M0 implementations.

## Non-Goals

- No FEM implementation.
- No SFL AST, diagnostics, regions, materials, elements, loads, or result IR.
- No mass rewrite of legacy headers.
- No ABI stability claim.
