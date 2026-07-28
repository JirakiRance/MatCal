# Compatibility Policy

MatCal M0 establishes a compatibility baseline. It does not declare production readiness or ABI stability.

## Compatibility Levels

Behavior compatibility means existing documented legacy behavior keeps working unless it is a confirmed bug with an explicit migration path.

Source/API compatibility means existing public headers, namespaces, class names, function names, constructor forms, and argument order remain available in the current major line.

Binary/ABI compatibility means already compiled consumers can link against a new binary without recompilation. MatCal does not promise this today. The current library is effectively header-only, compiler-specific, and not ABI-governed.

## M0 Commitments

- Do not delete public headers.
- Do not rename public namespaces or classes.
- Do not remove old functions directly.
- Prefer adding safe new APIs and keeping deprecated wrappers for public mistakes.
- Keep `MatCal::Legacy` as the target for the current public surface.
- Do not overwrite or replace historical compiled binaries used by PT or other projects.
- Do not promise cross-compiler or cross-version ABI stability.

## M0.1 Compatibility Notes

M0.1 is a source-compatible safety repair pass. It preserves public headers, namespaces, class names, and legacy entry points.

Source-compatible changes:

- `UpperTriangularMatrix` and `LowerTriangularMatrix` now use correct copy/move/assignment behavior while keeping the existing constructor forms available.
- `TridiagonalMatrix(0)` is now explicitly allowed; its internal bands are empty, and solving an empty right-hand side returns an empty result.
- `NumericalIntegration::Instant` now throws clear `std::invalid_argument` exceptions for an empty callable, non-finite endpoints, non-finite `eps`, and invalid `eps`.
- `OrthogonalPolynomials::Legendre` now returns mathematically correct recurrence coefficients.
- `determinant` now returns the correct sign when row swaps are required.
- `QinJiuShao::toFunction()` now captures polynomial coefficients by value, so the returned callable owns the state it needs.
- `SOR(AbstractMatrix&, AbstractMatrix&, double, ...)` is added as the precise overload. The old `int omega` overload remains and forwards to the double overload.
- `LU_Decompose` remains a no-pivot LU routine. Its failure message was clarified; pivoted LU is still a future API.

Behavior changes are intentional bug fixes. They may affect code that accidentally depended on broken behavior, such as `Instant` always clearing the callable, `determinant` returning the wrong sign after a row swap, or `toFunction()` reflecting mutations to a still-live polynomial after the callable was created.

ABI compatibility is not promised. Because MatCal remains header-only, existing consumers should recompile to pick up M0.1 fixes. M0.1 does not overwrite or replace any historical compiled MatCal binary used by PT or other projects.

## M1 Parallel API Policy

M1 adds `MatCal::Linalg` as a new 0.x development target while preserving `MatCal::Legacy`.

- `MatCal::Legacy` remains available and does not depend on `MatCal::Linalg`.
- `MatCal::Linalg` does not depend on legacy `AbstractMatrix`, `MatCal::Utils::Matrix`, dynamic casts, or legacy unique-pointer result patterns.
- No legacy header, namespace, class, or function is removed.
- Existing PT-style consumers can continue to use the legacy headers.
- New linalg consumers should include headers from `include/MatCal/Linalg` and link `MatCal::Linalg`.
- `MatCal::Linalg` is not ABI-stable in M1 and should be treated as a 0.x source API.

## M2/M2.1 Linalg Freezing Notes

M2 adds `SymmetricSkylineMatrix`, `SkylineLdltFactorization`, `factorize_skyline_ldlt`, and `solve_skyline_ldlt` under `MatCal::Linalg`. M2.1 keeps those existing API names and implementation path, then optimizes profile access and clarifies solver metrics.

- `MatCal::Legacy` remains unchanged and source-compatible with the M0/M0.1 baseline.
- `MatCal::Linalg` remains a 0.x source API and is not ABI-stable.
- M2.1 extends `SolverMetrics` with `factorization_operation_count` and `solve_operation_count`; existing `operation_count` remains available.
- M2.1 changes the 0.x tolerance composition to `max(absolute_tolerance, relative_tolerance * scale)`. This is a pre-1.0 contract correction for integration readiness.
- No precompiled MatCal binary is declared as the default integration vehicle. Consumers should build and link the CMake target for the pinned source version.

## M3 Legacy Migration Policy

M3 starts moving correct Legacy algorithms into shared modern cores.

- `MatCal::Legacy` remains available with the old headers, namespaces, class names, and main signatures.
- Legacy may depend on a reviewed shared core when doing so preserves source compatibility.
- `QinJiuShao` now delegates core polynomial algorithms to `MatCal::Polynomial::Polynomial`.
- Legacy output formatting and old node-list access stay in the facade when they are compatibility behavior rather than reusable math core.
- Installed packages include legacy headers at the include root so old includes such as `#include "QinJiuShao.hpp"` continue to work for target-based consumers.
- ABI stability is still not promised.

## M4 Legacy Migration Policy

M4 continues the facade strategy for scalar roots and interpolation.

- `Bisection`, `Picard`, `Newton`, and `Secant` keep their legacy class names, function names, argument order, old result structs, and exception-oriented `solve` behavior.
- `LinearInsert` and `CubicSpline` keep their legacy constructors, `calculate`, `getXs`, `getYs`, and `getM` accessors.
- Legacy scalar root APIs now delegate to `MatCal::Roots`.
- Legacy `LinearInsert` and `CubicSpline` now delegate to `MatCal::Interpolation`.
- The new interpolation core defaults to rejecting extrapolation, but legacy interpolation explicitly uses extrapolation to preserve old PT-style behavior.
- Bisection no longer treats a tiny absolute residual alone as endpoint success. This fixes scaled-function false positives such as `1e-20 * (x - c)`.
- Consumers using the new source package should recompile. ABI stability and compatibility with old precompiled binaries are still not promised.

## Deprecation Strategy

When an existing public API is unsafe but likely used:

1. Add a replacement with clearer ownership, constness, status reporting, and numeric contract.
2. Keep the legacy entry point forwarding to the replacement where feasible.
3. Mark the legacy API deprecated only after tests and migration notes exist.
4. Remove only at a declared major-version boundary.

## Version Boundary

Suggested boundaries:

- `0.x`: compatibility baseline, tests, internal refactor, no ABI promise.
- `1.x`: formal source-compatible stable legacy facade plus selected new APIs.
- `2.0`: first possible removal window for deprecated APIs, with migration guide.
