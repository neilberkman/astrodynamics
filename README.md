# astrodynamics

`astrodynamics` is a pure Rust numerical astrodynamics library for orbit
propagation, force models, and future flight-dynamics tooling.

The current crate focuses on a clean propagation foundation:

- Cartesian inertial state representation
- Two-body gravity
- J2 perturbation
- Fixed-step RK4
- Adaptive Dormand-Prince 5(4) (`DP54`)
- Oracle tests for Kepler invariants and J2 secular RAAN drift

This is a deliberately scoped `0.5.0`. It is a propagation core, not yet a
complete flight-dynamics platform. Dense output, event detection, richer
contexts, estimation, maneuvers, and higher-fidelity force models are planned
future work.

## Units and conventions

- Position: kilometers
- Velocity: kilometers per second
- Time: seconds
- State frame: inertial Cartesian state

## Example

```rust
use astrodynamics::forces::{CompositeForceModel, J2Gravity, TwoBodyGravity};
use astrodynamics::integrators::{DP54, Integrator};
use astrodynamics::propagator::{OrbitalDynamics, PropagationContext, api::IntegratorOptions};
use astrodynamics::CartesianState;

let initial = CartesianState::new(
    0.0,
    [7000.0, 0.0, 0.0],
    [0.0, 7.54605329, 0.0],
);

let mut forces = CompositeForceModel::new();
forces.add(Box::new(TwoBodyGravity::default()));
forces.add(Box::new(J2Gravity::default()));

let dynamics = OrbitalDynamics { force_model: &forces };
let integrator = DP54;
let ctx = PropagationContext::default();
let opts = IntegratorOptions {
    abs_tol: 1.0e-12,
    rel_tol: 1.0e-12,
    ..IntegratorOptions::default()
};

let result = integrator
    .propagate(initial, 3600.0, &dynamics, &ctx, &opts)
    .expect("propagation failed");

let final_state = result.final_state;
assert!(final_state.position_km.norm() > 6000.0);
```

## Validation

The current test suite includes:

- circular-orbit full-period return checks
- elliptic-orbit energy and angular-momentum invariants
- analytical J2 secular RAAN drift comparison

## Status

This crate is currently a focused propagation foundation, not yet a complete
mission-analysis suite. It provides inertial Cartesian state propagation with
RK4 and adaptive Dormand-Prince 5(4), plus two-body and J2 force models, with
additional flight-dynamics capabilities planned over time.

## License

MIT
