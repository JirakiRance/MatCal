# Warning Baseline

This document records warning status after M0.1. It is a compatibility tool, not a waiver for future APIs.

## Fixed in M0.1

- Triangular copy constructors no longer trigger base-class-copy initialization warnings in the regression path.
- `NumericalIntegration::Instant` no longer uses assignment as a truth value when checking the callable.

## Still Present in Legacy Headers

- Some loops compare signed `int` indices with unsigned `size_t` container sizes.
- `solve_Linear_System` is an inline legacy helper that may appear unused in focused test builds.

## M1 Handling

M1 should reduce these warnings through compatibility-preserving overloads, local type cleanup, and internal helper extraction. Do not batch-reformat or rename public APIs only to silence warnings.

## Compatibility Rule

Warnings that expose confirmed unsafe behavior should get tests and targeted fixes. Warnings tied to legacy shape, style, or naming should remain documented until a planned compatibility facade exists.

## M1 Linalg Warning Rule

`matcal_linalg` and its tests are built with stricter warning settings than legacy tests. On GCC/Clang-family compilers this includes `-Wall -Wextra -Wpedantic -Werror`.

Legacy warning debt remains documented here and is not required to be zeroed before M1 is accepted.

## Fixed in M4

- `QinJiuShao` now declares copy/move assignment explicitly. This removes the deprecated implicit copy-assignment warning source without changing the public class name or constructor set.
- The migrated `Picard::solve` path no longer owns the old unused `last_delta` local.

## Fixed in M5

- `CompositeNewtonCotes` no longer owns the old independent loop with an unused `seg_b` local; it delegates to `MatCal::Calculus`.

## M5 Remaining Legacy Warnings

- Multivariable derivative helpers and least-squares paths in `Basics.hpp` still compare old `int` dimensions with container `size()`. They remain legacy-only after M5.

## Fixed in M6

- `Insert.hpp` no longer returns scalar `getDegree()` values as `const int`.
- Legacy `NewtonForEquations` and `Least_Square` no longer own the migrated loops that contributed extra signed/unsigned warnings.

## M6 Remaining Legacy Warnings

- `QinJiuShao(int, std::vector<double>&)` still compares an `int` loop index with `zeros.size()`.
- `solve_Linear_System` may still appear unused in focused test builds.

## Fixed in M7

- `Derivative::pF_px` and `Derivative::dF_dx` no longer own their old finite-difference loops. They delegate to `MatCal::Calculus`, which also removes the previous `Derivative::dF_dx` signed/unsigned loop warning source.

## Fixed in M7.1

- `QinJiuShao(int, std::vector<double>&)` now uses an unsigned loop index when walking `zeros`.
- Deprecated private `sortByDegree()` is no longer called from the deprecated private `remove_fake()` helper body.
- `solve_Linear_System` is marked `[[maybe_unused]]` to preserve the legacy helper while avoiding focused-build unused warnings.

## M7.1 Remaining Legacy Warnings

- No known GCC warning remains from the M7.1 targeted legacy warning list.
