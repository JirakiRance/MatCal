# Getting Started

MatCal is a C++20 CMake project. For new code, link the smallest module target
you need instead of linking every target.

## Build and Test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## add_subdirectory

```cmake
add_subdirectory(MatCal)
target_link_libraries(app PRIVATE MatCal::Linalg)
```

## install and find_package

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
cmake --install build --prefix /path/to/matcal-install
```

```cmake
find_package(MatCal 0.3 CONFIG REQUIRED)
target_link_libraries(app PRIVATE MatCal::Linalg)
```

## Minimal Solver

```cpp
#include "MatCal/Linalg/DenseMatrix.hpp"
#include "MatCal/Linalg/DenseSolver.hpp"
#include "MatCal/Linalg/Vector.hpp"

#include <stdexcept>

int main() {
    MatCal::Linalg::DenseMatrix a{{4.0, 1.0}, {1.0, 3.0}};
    MatCal::Linalg::Vector b{1.0, 2.0};
    auto result = MatCal::Linalg::solve_dense_partial_pivot(a, b);
    if (!result.success()) {
        throw std::runtime_error("linear solve failed");
    }
}
```

This pattern is compiled by `matcal_documentation_example_smoke_tests`.
