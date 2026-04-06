//! Numerical astrodynamics engine for orbit propagation, force models, and
//! future flight-dynamics primitives.
//!
//! Current scope:
//!
//! - inertial Cartesian state representation
//! - two-body gravity and J2 perturbation
//! - fixed-step RK4
//! - adaptive Dormand-Prince 5(4) (`DP54`)
//! - propagation results with step statistics
//!
//! Planned future work includes dense output, event handling, richer propagation
//! contexts, additional force models, covariance propagation, estimation, and
//! maneuver support.

pub mod error;
pub mod state;
pub mod time;
pub mod constants;
pub mod math;
pub mod propagator;
pub mod integrators;
pub mod forces;
pub mod events;

pub use error::PropagationError;
pub use state::CartesianState;
pub use time::Time;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forces::TwoBodyGravity;
    use crate::integrators::{DP54, Integrator};
    use crate::propagator::{OrbitalDynamics, PropagationContext, api::IntegratorOptions};
    use nalgebra::Vector3;

    #[test]
    fn test_two_body_dp54_precision() {
        let r_mag: f64 = 7000.0;
        let mu: f64 = 398600.4418;
        let v_mag: f64 = (mu / r_mag).sqrt();
        let initial_state = CartesianState {
            epoch_tdb_seconds: 0.0,
            position_km: Vector3::new(r_mag, 0.0, 0.0),
            velocity_km_s: Vector3::new(0.0, v_mag, 0.0),
        };
        
        let force = TwoBodyGravity::default();
        let dynamics = OrbitalDynamics { force_model: &force };
        let integrator = DP54;
        let ctx = PropagationContext::default();
        let opts = IntegratorOptions {
            abs_tol: 1e-12,
            rel_tol: 1e-12,
            initial_step: 1.0,
            min_step: 1e-15,
            ..IntegratorOptions::default()
        };
        
        let period = 2.0 * std::f64::consts::PI * (r_mag.powi(3) / mu).sqrt();
        let result = integrator.propagate(initial_state, period, &dynamics, &ctx, &opts).unwrap();
        
        let final_pos = result.final_state.position_km;
        let final_vel = result.final_state.velocity_km_s;
        
        // Oracle 1: Return to start precision (Sub-millimeter)
        assert!((final_pos.x - r_mag).abs() < 1e-7, "Position X error too large: {}", (final_pos.x - r_mag).abs());
        assert!(final_pos.y.abs() < 1e-7, "Position Y error too large: {}", final_pos.y.abs());
        
        // Oracle 2: Energy conservation (Specific mechanical energy)
        let initial_energy = v_mag.powi(2) / 2.0 - mu / r_mag;
        let final_v_mag = final_vel.norm();
        let final_r_mag = final_pos.norm();
        let final_energy = final_v_mag.powi(2) / 2.0 - mu / final_r_mag;
        assert!((final_energy - initial_energy).abs() < 1e-10, "Energy conservation failure: {}", (final_energy - initial_energy).abs());
    }
}
