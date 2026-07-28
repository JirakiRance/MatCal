# CMake Package

M3 added installable CMake package metadata for MatCal. M4 extends the package exports to Roots and Interpolation. M5 extends the package exports to Calculus and ODE. M6 extends the package exports to Nonlinear and LeastSquares.

## Build-Tree Use

Existing consumers can continue to use:

```cmake
add_subdirectory(path/to/MatCal)
target_link_libraries(app PRIVATE MatCal::Linalg)
```

M3 also exposes:

```cmake
target_link_libraries(app PRIVATE MatCal::Polynomial)
target_link_libraries(app PRIVATE MatCal::Roots)
target_link_libraries(app PRIVATE MatCal::Interpolation)
target_link_libraries(app PRIVATE MatCal::Calculus)
target_link_libraries(app PRIVATE MatCal::ODE)
target_link_libraries(app PRIVATE MatCal::Nonlinear)
target_link_libraries(app PRIVATE MatCal::LeastSquares)
target_link_libraries(app PRIVATE MatCal::Legacy)
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

## Compatibility Notes

- The package does not promise ABI stability.
- The default integration model is source build through CMake, not a precompiled binary.
- SFL/FEM integration is still out of scope for M6.
