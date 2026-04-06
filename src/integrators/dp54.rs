use crate::state::{CartesianState, StateDerivative};
use crate::propagator::api::{PropagationContext, IntegratorOptions};
use crate::propagator::result::{PropagationResult, PropagationPoint, PropagationStats};
use crate::propagator::controller::PIController;
use crate::error::PropagationError;
use crate::integrators::{Integrator, DynamicsModel};
use crate::integrators::tableau::DP54Tableau;
use nalgebra::Vector3;

pub struct DP54;

impl Integrator for DP54 {
    fn propagate(
        &self,
        initial: CartesianState,
        t_end_seconds: f64,
        rhs: &dyn DynamicsModel,
        ctx: &PropagationContext,
        opts: &IntegratorOptions,
    ) -> Result<PropagationResult, PropagationError> {
        let tableau = DP54Tableau::default();
        let mut controller = PIController::default();
        controller.order = 5.0;
        
        let mut state = initial;
        let mut t = initial.epoch_tdb_seconds;
        let dt_target = t_end_seconds - t;
        let sign = dt_target.signum();
        let target_abs = dt_target.abs();
        
        let mut h = opts.initial_step.min(target_abs) * sign;
        let mut steps_accepted = 0;
        let mut steps_rejected = 0;
        let mut evals = 0;
        let mut points = Vec::new();

        points.push(PropagationPoint {
            epoch_tdb_seconds: t,
            position_km: state.position_array(),
            velocity_km_s: state.velocity_array(),
        });

        // FSAL: k1
        let mut k1 = rhs.derivative(&state, ctx)?;
        evals += 1;

        while (t - initial.epoch_tdb_seconds).abs() < target_abs {
            if steps_accepted + steps_rejected >= opts.max_steps {
                return Err(PropagationError::MaxStepsExceeded);
            }

            let mut h_step = h;
            if (t + h_step - initial.epoch_tdb_seconds).abs() > target_abs {
                h_step = t_end_seconds - t;
            }

            // Step using DP54
            let (next_state, k_fsal, r_err, v_err, step_evals) = self.step(state, h_step, k1, rhs, ctx, &tableau)?;
            
            // Error estimation
            let r_scale = opts.abs_tol + state.position_km.norm().max(next_state.position_km.norm()) * opts.rel_tol;
            let v_scale = opts.abs_tol + state.velocity_km_s.norm().max(next_state.velocity_km_s.norm()) * opts.rel_tol;
            
            let err_r = r_err.norm() / r_scale;
            let err_v = v_err.norm() / v_scale;
            let err = err_r.max(err_v);

            if err <= 1.0 {
                // Accepted
                state = next_state;
                t += h_step;
                k1 = k_fsal; // FSAL
                steps_accepted += 1;
                evals += step_evals;
                
                if opts.dense_output {
                    points.push(PropagationPoint {
                        epoch_tdb_seconds: t,
                        position_km: state.position_array(),
                        velocity_km_s: state.velocity_array(),
                    });
                }
                
                h = controller.next_step(h_step, err);
            } else {
                steps_rejected += 1;
                evals += step_evals - 1;
                h = controller.next_step(h_step, err);
                
                if h.abs() < opts.min_step {
                    return Err(PropagationError::NumericalFailure("Step size too small".to_string()));
                }
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
                accepted_steps: steps_accepted,
                rejected_steps: steps_rejected,
                evaluations: evals,
            },
        })
    }
}

impl DP54 {
    fn step(
        &self,
        state: CartesianState,
        h: f64,
        k1: StateDerivative,
        rhs: &dyn DynamicsModel,
        ctx: &PropagationContext,
        tableau: &DP54Tableau,
    ) -> Result<(CartesianState, StateDerivative, Vector3<f64>, Vector3<f64>, u32), PropagationError> {
        let mut ks = Vec::with_capacity(7);
        ks.push(k1);

        for i in 1..6 {
            let mut dpos = Vector3::zeros();
            let mut dvel = Vector3::zeros();
            for j in 0..i {
                dpos += ks[j].dpos_km_s * tableau.a[i][j];
                dvel += ks[j].dvel_km_s2 * tableau.a[i][j];
            }
            
            let stage_state = CartesianState {
                epoch_tdb_seconds: state.epoch_tdb_seconds + h * tableau.c[i],
                position_km: state.position_km + dpos * h,
                velocity_km_s: state.velocity_km_s + dvel * h,
            };
            ks.push(rhs.derivative(&stage_state, ctx)?);
        }

        // 5th order solution
        let mut dpos5 = Vector3::zeros();
        let mut dvel5 = Vector3::zeros();
        for i in 0..6 {
            dpos5 += ks[i].dpos_km_s * tableau.b5[i];
            dvel5 += ks[i].dvel_km_s2 * tableau.b5[i];
        }

        let next_state = CartesianState {
            epoch_tdb_seconds: state.epoch_tdb_seconds + h,
            position_km: state.position_km + dpos5 * h,
            velocity_km_s: state.velocity_km_s + dvel5 * h,
        };

        // FSAL
        let k_fsal = rhs.derivative(&next_state, ctx)?;
        ks.push(k_fsal);

        // 4th order for error estimate
        let mut dpos4 = Vector3::zeros();
        let mut dvel4 = Vector3::zeros();
        for i in 0..7 {
            dpos4 += ks[i].dpos_km_s * tableau.b4[i];
            dvel4 += ks[i].dvel_km_s2 * tableau.b4[i];
        }

        let r_err = (dpos5 - dpos4) * h;
        let v_err = (dvel5 - dvel4) * h;

        Ok((next_state, k_fsal, r_err, v_err, 6))
    }
}
