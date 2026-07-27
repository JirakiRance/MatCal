# MatCal M0 Code Audit

Audit date: 2026-07-27.

## Repository State

- Branch: `main`.
- HEAD: `ce295a7a39cbbd21c1a39dde4daedf9180aa2a0d` (`添加了readme`).
- Remotes:
  - `github`: `https://github.com/JirakiRance/MatCal.git`
  - `origin`: `https://gitee.com/jirakirance/mat-cal.git`
- Initial status: clean relative to `github/main`. Git emitted permission warnings for `C:\Users\DELL/.config/git/ignore`; repository state itself was clean.
- Tracked files at audit start: `.gitignore`, `README.md`, `src/Basics.hpp`, `src/Insert.hpp`, `src/Iteration.hpp`, `src/Matrix.hpp`, `src/QinJiuShao.hpp`.
- Untracked files at audit start: none reported by `git status --short --branch`.
- EOL baseline: `git ls-files --eol` reported `i/lf w/crlf` for every tracked file. This means the index stores LF while the working tree uses CRLF. No pure line-ending changes should be committed.

## Recent History

Recent commits show MatCal is an existing handwritten numerical library, not a greenfield project:

- `ce295a7` add README.
- `fbfc8c2` ignore build/editor outputs and stop tracking them.
- `bc2cdb9` end-of-course snapshot and bug fixes.
- `d1b9fae` Matrix parent-pointer constructor; orthogonal polynomials, numerical integration, ODE.
- `1f1c1bb` PT project updated spline interpolation and linear interpolation; add integration, RK4, Euler integration.
- `a2877c7` move least squares into `Least_Square::solve`.
- `a734a50` implement `TridiagonalMatrix` and chase/Thomas solver.
- Earlier commits add nonlinear equations, interpolation namespaces, SOR, Jacobi, Gauss-Seidel, LU and QinJiuShao.

## README Findings

`README.md` is an API manual for the current legacy library. It documents:

- `MatCal::Utils`: `QinJiuShaoNode`, `QinJiuShao`, `AbstractMatrix`, `Matrix`, `SparseMatrix`, upper/lower triangular matrices, `TridiagonalMatrix`.
- `MatCal::Algorithm::Matrix`: elementary row/column transformations, norms, Gaussian elimination, LU, Jacobi, Gauss-Seidel, SOR, power methods.
- `MatCal::Algorithm::Basics`: derivatives, least squares, orthogonal polynomials, numerical integration, ODE.
- `MatCal::Algorithm::Insert`: Lagrange, Newton quotient, Newton finite-difference, Hermite, cubic spline, linear interpolation.
- `MatCal::Algorithm::Iteration`: scalar and vector nonlinear solvers.

The README does not explicitly describe SFL. The Git history and source comments explicitly mention PT, especially interpolation and `MatCal::Algorithm::Basics::Integrate::RK4`.

## Build Baseline Added

This M0 pass added standard CMake support:

- `MatCal::Legacy`: current header-only legacy API target.
- C++20 compile features.
- CTest tests in `tests/`.
- GCC/Clang/MSVC warning options for test targets.
- Optional `MATCAL_ENABLE_SANITIZERS` for GCC/Clang-family builds.

No public header was removed or renamed.

## Actual Build/Test Result

Observed platform:

- OS shell: Windows PowerShell in `E:\Ducuments\codeForVScode\MatCal`.
- CMake: 4.4.0.
- Compiler used by CMake: GNU C++ 15.1.0, `E:/minGW/mingw64/bin/c++.exe`.
- Clang/MSVC: not found in PATH during this audit.

Commands run:

```powershell
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Debug
cmake --build build --config Debug
ctest --test-dir build --output-on-failure
```

CTest result: 2/2 tests passed.

Warnings remain. They are documented rather than batch-fixed in M0 because many are tied to legacy copy rules, signedness, and numerical contracts.

## Engineering Changes Made

- Added `CMakeLists.txt`.
- Added characterization tests.
- Added multi-translation-unit link test.
- Marked header-defined non-template functions/static data as `inline` where required to make the legacy header-only library link safely from multiple translation units.
- Added missing standard includes needed by the existing header definitions.

No algorithmic rewrite was performed.
