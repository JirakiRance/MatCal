# MatCal Tests

`tests/` is the tracked, reproducible CTest suite used by local builds and CI.
Do not add this directory to `.gitignore`.

## Layout

- `unit/legacy/`: characterization and compatibility tests for the legacy facade.
- `unit/linalg/`: Linalg container, solver, Skyline/LDLT, pivoted LU, iterative, and eigen tests.
- `unit/<module>/`: focused tests for Polynomial, Roots, Interpolation, Calculus, ODE, Nonlinear, and LeastSquares.
- `integration/public_headers/`: public-header self-contained compile tests.
- `integration/multi_tu/`: multi-translation-unit ODR/link tests.
- `integration/consumers/`: add_subdirectory, install/find_package, SFL-compatible, and PT-style consumers.
- `integration/docs_smoke/`: examples mirrored from user documentation.

The ignored `test/` directory is separate. It contains historical handwritten
personal experiments and is not part of CI. Generated binaries under
`test/output/` may be deleted because they are reproducible build products.
