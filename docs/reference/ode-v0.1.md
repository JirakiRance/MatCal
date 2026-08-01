# MatCal::ODE v0.1

M5 adds `MatCal::ODE` as an independent header-only target for value-safe ordinary differential equation stepping. It is domain-neutral and does not contain PT, SFL, FEM, force, material, mesh, or result-recovery concepts.

## Public Header

```cpp
#include <MatCal/ODE/ODE.hpp>
```

## State Contract

The new core represents state as `std::vector<double>` and RHS callables as:

```cpp
std::function<std::vector<double>(double t, const std::vector<double>& y)>
```

All state values, time values, step sizes, and RHS outputs must be finite. RHS output size must match the input state size. Empty states and size mismatches are rejected with structured diagnostics.

## Result Contract

`OdeStepResult` owns:

- output state;
- output time;
- completion flag;
- machine-readable `OdeDiagnostic`;
- RHS evaluation count.

`OdeTrajectoryResult` owns the trajectory rows. Rows use the legacy-compatible shape `[t, y0, y1, ...]`.

Numerical failure does not return a pseudo-success trajectory.

## RK4

`rk4_step` reuses the legacy classic four-stage RK4 formula from `ODE::RungeKutta_44` and `Integrate::RK4::step`:

```text
k1 = f(t, y)
k2 = f(t + h/2, y + h*k1/2)
k3 = f(t + h/2, y + h*k2/2)
k4 = f(t + h,   y + h*k3)
y_next = y + h*(k1 + 2*k2 + 2*k3 + k4)/6
```

Negative `dt` is allowed by default for the new core and is documented as backward stepping. Legacy ODE table APIs still require positive `h` to preserve their old source contract.

## Euler

`euler_step` reuses the legacy simple Euler formula. `improved_euler_step` reuses the legacy predictor/corrector formula from `ODE::Euler`.

## Legacy Delegation

The following legacy APIs now forward to this core:

- `ODE::SimpleEuler`;
- `ODE::Euler`;
- `ODE::RungeKutta_44`;
- `Integrate::RK4::step`;
- `Integrate::RK4::step2`.

The PT-sensitive `RK4::step` and `step2` signatures remain available. Existing precompiled binaries are not replaced; source consumers should recompile to pick up M5 behavior.
