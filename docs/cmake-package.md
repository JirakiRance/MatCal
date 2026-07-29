# CMake Package

M3 added installable CMake package metadata for MatCal. M4 extends the package exports to Roots and Interpolation. M5 extends the package exports to Calculus and ODE. M6 extends the package exports to Nonlinear and LeastSquares. M7 adds more public headers under the existing `MatCal::Linalg` and `MatCal::Calculus` targets. M8 adds reusable pivoted LU under `MatCal::Linalg`. M8.1 prepares the package metadata for the `v0.3.0-alpha.1` release candidate.

The CMake project version is `0.3.0`. The release candidate tag is expected to be named `v0.3.0-alpha.1`; CMake package version files use the numeric `0.3.0` value. MatCal remains a 0.x source API and does not promise ABI stability.

## Build-Tree Use

Existing consumers can continue to use:

```cmake
add_subdirectory(path/to/MatCal)
target_link_libraries(app PRIVATE MatCal::Linalg)
```

Current exported targets are:

```cmake
target_link_libraries(app PRIVATE MatCal::Legacy)
target_link_libraries(app PRIVATE MatCal::Linalg)
target_link_libraries(app PRIVATE MatCal::Polynomial)
target_link_libraries(app PRIVATE MatCal::Roots)
target_link_libraries(app PRIVATE MatCal::Interpolation)
target_link_libraries(app PRIVATE MatCal::Calculus)
target_link_libraries(app PRIVATE MatCal::ODE)
target_link_libraries(app PRIVATE MatCal::Nonlinear)
target_link_libraries(app PRIVATE MatCal::LeastSquares)
```

M7 iterative solvers and power eigen solvers are used through `MatCal::Linalg`:

```cpp
#include "MatCal/Linalg/IterativeSolvers.hpp"
#include "MatCal/Linalg/EigenSolvers.hpp"
```

M7 multivariable derivative helpers are used through `MatCal::Calculus`:

```cpp
#include "MatCal/Calculus/Calculus.hpp"
```

## Install-Tree Use

After installing MatCal, consumers can use:

```cmake
find_package(MatCal CONFIG REQUIRED)
target_link_libraries(app PRIVATE
    MatCal::Legacy
    MatCal::Linalg
    MatCal::Polynomial
    MatCal::Roots
    MatCal::Interpolation
    MatCal::Calculus
    MatCal::ODE
    MatCal::Nonlinear
    MatCal::LeastSquares
)
```

Installed files include:

- `MatCalTargets.cmake`
- `MatCalConfig.cmake`
- `MatCalConfigVersion.cmake`
- public headers under `include/MatCal/...`
- legacy headers at the include root for source compatibility with existing includes such as `#include "QinJiuShao.hpp"`

## Package Review

M8.1 package consumers check that all formal targets are exported:

- `MatCal::Legacy`
- `MatCal::Linalg`
- `MatCal::Polynomial`
- `MatCal::Roots`
- `MatCal::Interpolation`
- `MatCal::Calculus`
- `MatCal::ODE`
- `MatCal::Nonlinear`
- `MatCal::LeastSquares`

The package consumer links only the targets; it does not add manual include directories. This checks that installed include paths and transitive dependencies propagate through targets. `MatCal::Legacy` exports its dependencies on migrated header-only cores. `MatCal::Linalg` remains independent from `MatCal::Legacy`.

The generated `MatCalConfigVersion.cmake` uses `SameMinorVersion` compatibility. This is intentional for the 0.x line: a consumer requesting `0.3` may accept compatible `0.3.z` packages, while `0.2` and `0.4` are not treated as interchangeable.

## Compatibility Notes

- The package does not promise ABI stability.
- The default integration model is source build through CMake, not a precompiled binary.
- SFL/FEM implementation is still out of scope for M8.1.
- SFL may pin and build this package through CMake during a controlled upgrade from `v0.2.0-alpha.1` to `v0.3.0-alpha.1`.
