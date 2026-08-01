# MatCal

[![CI](https://github.com/JirakiRance/MatCal/actions/workflows/ci.yml/badge.svg)](https://github.com/JirakiRance/MatCal/actions/workflows/ci.yml)

MatCal is a standalone C++20 numerical calculation library. It began as a
handwritten legacy C++ library and now provides modern CMake targets while
keeping the original headers and namespaces available for older projects.

Current release: `v0.3.0-alpha.2`

MatCal is in maintenance mode for the current project stage. The completed
scope is enough for the present SFL/course-design use cases, especially dense
linear algebra, symmetric Skyline SPD LDLT, pivoted LU, and structured solver
results. New algorithm work is paused until a real downstream need or confirmed
bug appears.

## Modules

- `MatCal::Linalg`: `Vector`, `DenseMatrix`, partial-pivot LU, Skyline SPD LDLT,
  stationary iterative solvers, and power-method eigen solvers.
- `MatCal::Polynomial`: owning dense polynomial utilities.
- `MatCal::Roots`: scalar root solvers with structured results.
- `MatCal::Interpolation`: linear, cubic spline, Lagrange, Newton, and Hermite
  interpolation.
- `MatCal::Calculus`: finite differences and numerical integration.
- `MatCal::ODE`: Euler and RK4 helpers.
- `MatCal::Nonlinear`: multivariable Newton solver.
- `MatCal::LeastSquares`: polynomial least-squares fitting.
- `MatCal::Legacy`: compatibility facade for historical headers under `src/`.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## add_subdirectory

```cmake
add_subdirectory(MatCal)
target_link_libraries(your_target PRIVATE MatCal::Linalg)
```

## install and find_package

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
cmake --install build --prefix /path/to/matcal-install
```

```cmake
find_package(MatCal 0.3 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE MatCal::Linalg)
```

All exported targets are documented in
[docs/index.md](docs/index.md).

## Minimal Linalg Example

```cpp
#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/Vector.hpp"

#include <stdexcept>

int main() {
    MatCal::Linalg::DenseMatrix a{{4.0, 1.0}, {1.0, 3.0}};
    MatCal::Linalg::Vector b{1.0, 2.0};

    MatCal::Linalg::SolverResult result =
        MatCal::Linalg::solve_dense_partial_pivot(a, b);

    if (!result.success()) {
        throw std::runtime_error(result.diagnostics.empty()
            ? "MatCal solve failed"
            : result.diagnostics.front().message);
    }

    return result.solution.size() == 2 ? 0 : 1;
}
```

The example style above is compiled by the documentation smoke test in
`tests/integration/docs_smoke/`.

## SolverResult Handling

Modern solvers return structured result objects. Check `success()` or the
status enum before using the solution. Numeric failures such as singular
matrices, non-positive-definite Skyline systems, non-finite intermediate values,
and non-convergence are reported through status, diagnostics, and metrics rather
than by printing to stdout.

See [docs/reference/solver-result-contract.md](docs/reference/solver-result-contract.md)
and [docs/reference/numerical-contracts.md](docs/reference/numerical-contracts.md).

## Legacy Compatibility

Legacy headers remain installed:

- `Basics.hpp`
- `Insert.hpp`
- `Iteration.hpp`
- `Matrix.hpp`
- `QinJiuShao.hpp`

They exist because PT and older experiments may still include historical
headers, namespaces, class names, and constructor forms. MatCal does not remove
Legacy simply to make the tree look modern; compatible facades are the migration
strategy. New code should prefer the module targets above.

No cross-version or cross-compiler ABI stability is promised for the alpha
series. Pin a Git tag and rebuild from source.

## SFL Boundary

SFL depends on MatCal, not the other way around:

```text
SFL -> MatCal
```

MatCal must not contain SFL AST, CAE IR, mesh/region, FEM element, material,
load, diagnostic, or result-recovery semantics. SFL should link the required
CMake target, currently primarily `MatCal::Linalg`, and should not copy MatCal
source files into its own repository.

## Documentation

Start at [docs/index.md](docs/index.md).

Important entries:

- [Getting Started](docs/getting-started.md)
- [Linalg User Guide](docs/user-guide/linalg.md)
- [Skyline LDLT User Guide](docs/user-guide/skyline-ldlt.md)
- [Pivoted LU User Guide](docs/user-guide/pivoted-lu.md)
- [Legacy Migration](docs/migration/legacy-migration.md)
- [Known Limitations](docs/reference/known-limitations.md)

## Current Limits

- `v0.3.0-alpha.2` is an alpha source API snapshot, not an ABI-stable release.
- Skyline LDLT is only for real symmetric positive definite matrices.
- Pivoted LU is dense partial pivoting, not QR or SVD.
- CSR, general sparse iterative solvers, QR, SVD, and indefinite symmetric
  factorization are intentionally out of scope for the current maintenance
  baseline.

## License

MatCal is licensed under the MIT License. See [LICENSE](LICENSE).
