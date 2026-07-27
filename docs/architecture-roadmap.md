# Architecture Roadmap

M0 keeps the legacy surface intact and makes future change measurable. M0.1 adds targeted safety repairs and regression tests without a full rewrite. M1 introduces the first independent `MatCal::Linalg` 0.x API.

## Current Shape

MatCal now has two parallel surfaces:

- `MatCal::Legacy`: the header-only legacy library with broad public APIs in five headers.
- `MatCal::Linalg`: the new 0.x value-type linalg target under `include/MatCal/Linalg`.

The two targets do not depend on each other in M1.

## Direction

Target layering:

```text
legacy MatCal API
    -> compatibility facade
    -> internal linalg/numerics implementation
```

Public targets:

- `MatCal::Legacy`: current API surface.
- `MatCal::Linalg`: 0.x development API for independent vector/matrix storage and reference solvers.

Future public targets may include:

- `MatCal::Core`: scalar utilities and shared tolerance/status policies if they become useful outside linalg.

Do not describe `MatCal::Linalg` as ABI-stable or production-ready.

## Gradual Refactor Plan

1. Keep M0/M0.1 characterization and regression tests as the compatibility gate.
2. Maintain M1 `Vector`, `DenseMatrix`, `SolverOptions`, and `SolverResult` as the new baseline.
3. Add pivoted LU factors as a new explicit API while preserving legacy no-pivot `LU_Decompose`.
4. Separate matrix storage from solver algorithms.
5. Replace stdout from core paths with structured status or caller-provided diagnostics.
6. Expand scale-aware tolerance helpers only when algorithms need them.
7. Forward legacy APIs to new internals one family at a time, behind compatibility tests.
8. Add CSR and iterative solver infrastructure after dense behavior is stable.
9. Add `SymmetricSkylineMatrix` and `LDLT` in M2 or later for future structural mechanics needs, without adding FEM semantics.

## M1 Baseline

M1 builds on the M0.1 safety baseline rather than reopening the legacy headers wholesale:

- `Vector` and row-major `DenseMatrix` exist as independent value types.
- `SolverStatus`, `SolverOptions`, `SolverDiagnostic`, `SolverMetrics`, and `SolverResult` exist.
- `solve_dense_partial_pivot` is a reference dense solver with structured failure.
- Public linalg headers are tested for self-contained includes and multi-TU use.
- Legacy APIs are still separate and unchanged by M1.

M2 should add Skyline/LDLT only after the M1 value, tolerance, and result contracts are stable enough to support them.

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
- No CSR, Skyline, LDLT, or FEM implementation in M1.
