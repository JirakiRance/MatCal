# Architecture Roadmap

M0 keeps the legacy surface intact and makes future change measurable. M0.1 adds targeted safety repairs and regression tests without a full rewrite.

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

1. Keep M0/M0.1 characterization and regression tests as the compatibility gate.
2. Introduce `Vector`, `DenseMatrix`, and `SolverStatus/SolverResult` internally.
3. Add pivoted LU as a new explicit API while preserving legacy no-pivot `LU_Decompose`.
4. Separate matrix storage from solver algorithms.
5. Replace stdout from core paths with structured status or caller-provided diagnostics.
6. Add scale-aware tolerance helpers.
7. Forward legacy APIs to new internals one family at a time.
8. Add CSR and iterative solver infrastructure after dense behavior is stable.
9. Add `SymmetricSkylineMatrix` and `LDLT` for future structural mechanics needs, without adding FEM semantics.

## M1 Starting Point

M1 should build on the M0.1 safety baseline rather than reopening the legacy headers wholesale:

- Add status/result types for direct and iterative solvers.
- Add const-safe overloads for copy-like matrix APIs.
- Add checked dense conversion helpers to remove scattered `dynamic_cast` assumptions.
- Add pivoted LU under a new name.
- Start forwarding one legacy solver family through an internal implementation while preserving legacy signatures.
- Keep `MatCal::Linalg` unpublished until storage, tolerance, and solver contracts are documented and tested.

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
