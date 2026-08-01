# MatCal::Polynomial v0.1

M3 adds `MatCal::Polynomial` as the first non-linalg modernized core extracted from Legacy behavior.

## CMake Target

- Library target: `matcal_polynomial`
- Alias target: `MatCal::Polynomial`

The target is header-only in M3. `MatCal::Legacy` links it as an interface dependency so installed legacy headers can include the shared polynomial core.

## Public Header

- `include/MatCal/Polynomial/Polynomial.hpp`

The header is self-contained and covered by a multi-translation-unit test.

## Value Type

`MatCal::Polynomial::Polynomial` is an owning value type. It stores coefficients densely in ascending degree order:

```text
coefficients()[i] = coefficient of x^i
```

This representation was chosen for the first migrated core because Horner evaluation, derivative, integral, and multiplication are simple and directly verifiable. Legacy `QinJiuShao` still exposes its sparse descending node facade and adapts to/from this core. High sparse degrees are tested, but a sparse polynomial core remains a future design option if real workloads require it.

## Contract

- The zero polynomial stores no coefficients.
- `degree()` returns `0` for the zero polynomial.
- Coefficients with absolute value below `1e-12` are trimmed.
- Constructors reject non-finite coefficients.
- `evaluate()` rejects non-finite input and throws if evaluation overflows to NaN/Inf.
- Scalar division by near-zero throws `std::invalid_argument`.
- `to_function()` captures a polynomial value copy, so the callable owns its state.

## Operations

M3 provides:

- `from_terms`
- `coefficient`
- `degree`
- `coefficients`
- `terms_descending`
- `evaluate`
- `derivative`
- `integral`
- `definite_integral`
- `to_function`
- polynomial add/subtract/multiply
- scalar multiply/divide

## Legacy Delegation

`MatCal::Utils::QinJiuShao` now delegates the following algorithms to `MatCal::Polynomial::Polynomial`:

- `calculate`
- `toFunction`
- `operator+`
- `operator-`
- `operator*`
- scalar `operator*`
- scalar `operator/`
- `derivative`
- `integral`
- `definiteIntegral`

The old output methods and node-list facade remain legacy-owned for source compatibility.
