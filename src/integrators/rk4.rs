use crate::state::{CartesianState, StateDerivative};
use crate::propagator::api::{PropagationContext, IntegratorOptions};
use crate::propagator::result::{PropagationResult, PropagationPoint, PropagationStats};
use crate::error::PropagationError;
use crate::integrators::{Integrator, DynamicsModel};

pub struct RK4;

impl Integrator for RK4 {
    fn propagate(
        &self,
        initial: CartesianState,
        t_end_seconds: f64,
        rhs: &dyn DynamicsModel,
        ctx: &PropagationContext,
        opts: &IntegratorOptions,
    ) -> Result<PropagationResult, PropagationError> {
        let mut state = initial;
        let mut t = initial.epoch_tdb_seconds;
        let dt_target = t_end_seconds - t;
        let sign = dt_target.signum();
        let target_abs = dt_target.abs();
        
        let h_initial = opts.initial_step.min(target_abs) * sign;
        let mut h = h_initial;
        let mut steps = 0;
        let mut points = Vec::new();
        
        points.push(PropagationPoint {
            epoch_tdb_seconds: t,
            position_km: state.position_array(),
            velocity_km_s: state.velocity_array(),
        });

        while (t - initial.epoch_tdb_seconds).abs() < target_abs {
            if steps >= opts.max_steps {
                return Err(PropagationError::MaxStepsExceeded);
            }
            
            if (t + h - initial.epoch_tdb_seconds).abs() > target_abs {
                h = t_end_seconds - t;
            }

            let next_state = self.step(state, h, rhs, ctx)?;
            state = next_state;
            t += h;
            steps += 1;

            if opts.dense_output {
                points.push(PropagationPoint {
                    epoch_tdb_seconds: t,
                    position_km: state.position_array(),
                    velocity_km_s: state.velocity_array(),
                });
            }
        }

        if !opts.dense_output {
             points.push(PropagationPoint {
                epoch_tdb_seconds: t,
                position_km: state.position_array(),
                velocity_km_s: state.velocity_array(),
            });
        }

        Ok(PropagationResult {
            final_state: state,
            points,
            events: Vec::new(),
            stats: PropagationStats {
                accepted_steps: steps,
                rejected_steps: 0,
                evaluations: steps * 4,
            },
            dense: None,
        })
    }
}

impl RK4 {
    fn step(
        &self,
        state: CartesianState,
        h: f64,
        rhs: &dyn DynamicsModel,
        ctx: &PropagationContext,
    ) -> Result<CartesianState, PropagationError> {
        let k1 = rhs.derivative(&state, ctx)?;
        
        let s2 = self.advance(&state, &k1, h / 2.0);
        let k2 = rhs.derivative(&s2, ctx)?;
        
        let s3 = self.advance(&state, &k2, h / 2.0);
        let k3 = rhs.derivative(&s3, ctx)?;
        
        let s4 = self.advance(&state, &k3, h);
        let k4 = rhs.derivative(&s4, ctx)?;
        
        let dpos = (k1.dpos_km_s + k2.dpos_km_s * 2.0 + k3.dpos_km_s * 2.0 + k4.dpos_km_s) * (h / 6.0);
        let dvel = (k1.dvel_km_s2 + k2.dvel_km_s2 * 2.0 + k3.dvel_km_s2 * 2.0 + k4.dvel_km_s2) * (h / 6.0);
        
        Ok(CartesianState {
            epoch_tdb_seconds: state.epoch_tdb_seconds + h,
            position_km: state.position_km + dpos,
            velocity_km_s: state.velocity_km_s + dvel,
        })
    }

    fn advance(&self, state: &CartesianState, deriv: &StateDerivative, h: f64) -> CartesianState {
        CartesianState {
            epoch_tdb_seconds: state.epoch_tdb_seconds + h,
            position_km: state.position_km + deriv.dpos_km_s * h,
            velocity_km_s: state.velocity_km_s + deriv.dvel_km_s2 * h,
        }
    }
}
