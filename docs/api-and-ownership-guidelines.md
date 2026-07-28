# API and Ownership Guidelines

These guidelines apply to new APIs and gradual refactors. Legacy APIs remain available through `MatCal::Legacy`.

## Values, References, and Pointers

- Prefer value types for small scalar results and matrix/vector handles with RAII storage.
- Prefer `const T&` for read-only inputs.
- Use non-const `T&` only when the function intentionally mutates input.
- Avoid raw owning pointers.
- Avoid `double**` in new APIs. Keep existing raw-pointer constructors as compatibility wrappers that copy input.
- Use `std::unique_ptr` only when dynamic polymorphism is required.

## Result Reporting

New solvers should return structured results:

```cpp
enum class SolverStatus {
    success,
    invalid_input,
    dimension_mismatch,
    non_finite_input,
    singular,
    not_positive_definite,
    breakdown,
    not_converged
};

struct SolverResult {
    SolverStatus status;
    Vector solution;
    SolverMetrics metrics;
    std::vector<SolverDiagnostic> diagnostics;
};
```

Do not encode failure only through stdout or magic values such as `error = -1`.

M1 `MatCal::Linalg` solvers return `SolverResult` for ordinary numerical failure. M4/M5 scalar numerical cores follow the same direction with structured result objects. Bounds errors, ragged matrix construction, size overflow, invalid object construction, and legacy compatibility wrappers may throw exceptions.

## Const-Correctness

Legacy signatures frequently use non-const references even when the function copies inputs. New APIs should make mutation explicit. Compatibility wrappers can keep old signatures and forward internally.

## Lifetime

- Returned matrices must own their storage.
- Views must explicitly document borrowed lifetime.
- A returned callable must not capture `this` unless the API name and docs make the borrowed lifetime obvious.
- `Vector::span()` and `DenseMatrix::row()` return borrowed views. They are invalidated by object destruction and by operations that reallocate storage.

## M1 Linalg Ownership

- `Vector` and `DenseMatrix` are owning value types with default copy/move behavior.
- New APIs do not expose owning raw pointers.
- New APIs do not return `std::unique_ptr<AbstractMatrix>`.
- `MatCal::Linalg` does not inherit from legacy `AbstractMatrix`.

## M3 Polynomial Ownership

- `MatCal::Polynomial::Polynomial` owns its coefficient storage.
- `to_function()` captures a value copy and does not borrow `this`.
- Legacy `QinJiuShao` adapts to/from the polynomial value type when using migrated math operations.

## M4 Roots and Interpolation Ownership

- `MatCal::Roots` accepts `std::function<double(double)>` by value/reference at the call boundary and does not store callables after the solve returns.
- `RootResult` owns its diagnostic, metrics, and optional iteration series.
- `MatCal::Interpolation::LinearInterpolator` and `CubicSpline` own their node and coefficient storage.
- Legacy `LinearInsert` and `CubicSpline` keep old getters but mirror data owned by the new interpolation core.

## M5 Calculus and ODE Ownership

- `MatCal::Calculus` accepts scalar callables at the call boundary and does not retain them after evaluation.
- `DerivativeResult` and `IntegrationResult` own their diagnostics, metrics, and optional Romberg table data.
- `MatCal::ODE` accepts RHS callables at the call boundary and does not retain them after stepping or integration.
- ODE states and trajectories are value-owned `std::vector<double>` containers.
- Legacy ODE and PT-style RK4 wrappers adapt their old callable shapes to the new value-state core.

## Printing and Diagnostics

Core library code should not print to stdout in new APIs. Use return status, exceptions for programmer errors, or optional caller-provided diagnostics.

## Exceptions

Use exceptions for invalid dimensions, out-of-range access, and invalid options. Numerical failure should prefer `SolverResult` with a non-success status in new solver APIs.
