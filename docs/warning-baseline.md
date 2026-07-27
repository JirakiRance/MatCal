# Warning Baseline

This document records warning status after M0.1. It is a compatibility tool, not a waiver for future APIs.

## Fixed in M0.1

- Triangular copy constructors no longer trigger base-class-copy initialization warnings in the regression path.
- `NumericalIntegration::Instant` no longer uses assignment as a truth value when checking the callable.

## Still Present in Legacy Headers

- `QinJiuShao` has deprecated implicit copy-assignment warnings because it defines a copy constructor but not a matching assignment operator.
- Some loops compare signed `int` indices with unsigned `size_t` container sizes.
- `Insert.hpp` has `const` on a scalar return type.
- `CompositeNewtonCotes` keeps an unused local in one path.
- `solve_Linear_System` is an inline legacy helper that may appear unused in focused test builds.
- `Picard::solve` keeps an unused `last_delta` local.

## M1 Handling

M1 should reduce these warnings through compatibility-preserving overloads, local type cleanup, and internal helper extraction. Do not batch-reformat or rename public APIs only to silence warnings.

## Compatibility Rule

Warnings that expose confirmed unsafe behavior should get tests and targeted fixes. Warnings tied to legacy shape, style, or naming should remain documented until a planned compatibility facade exists.

## M1 Linalg Warning Rule

`matcal_linalg` and its tests are built with stricter warning settings than legacy tests. On GCC/Clang-family compilers this includes `-Wall -Wextra -Wpedantic -Werror`.

Legacy warning debt remains documented here and is not required to be zeroed before M1 is accepted.
