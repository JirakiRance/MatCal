# Architecture Roadmap

M0 keeps the legacy surface intact and makes future change measurable. M0.1 adds targeted safety repairs and regression tests without a full rewrite. M1 introduces the first independent `MatCal::Linalg` 0.x API. M1.1 hardens its scale, finite-value, and result contracts. M2 adds general symmetric skyline storage and SPD LDLT. M2.1 optimizes and freezes that existing skyline baseline for later integration work.

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
2. Maintain M1.1 `Vector`, `DenseMatrix`, `SolverOptions`, `SolverDiagnostic`, `SolverMetrics`, and `SolverResult` as the new baseline.
3. Maintain M2/M2.1 `SymmetricSkylineMatrix` and SPD `SkylineLdltFactorization` as general numerical facilities.
4. Add pivoted LU factors as a new explicit API while preserving legacy no-pivot `LU_Decompose`.
5. Separate matrix storage from solver algorithms.
6. Replace stdout from core paths with structured status or caller-provided diagnostics.
7. Expand scale-aware tolerance helpers only when algorithms need them.
8. Forward legacy APIs to new internals one family at a time, behind compatibility tests.
9. Add CSR and iterative solver infrastructure after dense and skyline behavior is stable.

## M1 Baseline

M1 builds on the M0.1 safety baseline rather than reopening the legacy headers wholesale:

- `Vector` and row-major `DenseMatrix` exist as independent value types.
- `SolverStatus`, `SolverOptions`, `SolverDiagnostic`, `SolverMetrics`, and `SolverResult` exist.
- `solve_dense_partial_pivot` is a reference dense solver with structured failure.
- Public linalg headers are tested for self-contained includes and multi-TU use.
- Legacy APIs are still separate and unchanged by M1.

M2 reuses:

- Scale-relative tolerance form.
- Explicit `matrix_scale`, `rhs_scale`, `solution_scale`, and residual metrics.
- `not_positive_definite` vs `breakdown` distinction for SPD LDLT.
- Machine-readable diagnostic `code`, `reason`, `phase`, `row`, and `column`.

M2.1 adds:

- Profile-limited skyline back substitution using factorization-owned column adjacency.
- Separate factorization and solve work counters in `SolverMetrics`.
- Performance contract tests based on deterministic storage/work counts, not wall-clock timing.
- The frozen 0.x tolerance composition `max(absolute_tolerance, relative_tolerance * scale)`.

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
- No CSR, Skyline, LDLT, or FEM implementation in M1.1.
- No general indefinite symmetric factorization, Bunch-Kaufman, CSR, iterative method, FEM, or SFL integration in M2/M2.1.
