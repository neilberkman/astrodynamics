use crate::state::{CartesianState, StateDerivative};
use crate::propagator::api::PropagationContext;
use crate::error::PropagationError;
use crate::integrators::DynamicsModel;
use crate::forces::r#trait::ForceModel;

pub struct OrbitalDynamics<'a> {
    pub force_model: &'a dyn ForceModel,
}

impl<'a> DynamicsModel for OrbitalDynamics<'a> {
    fn derivative(
        &self,
        state: &CartesianState,
        ctx: &PropagationContext,
    ) -> Result<StateDerivative, PropagationError> {
        let accel = self.force_model.acceleration(state, ctx)?;
        Ok(StateDerivative {
            dpos_km_s: state.velocity_km_s,
            dvel_km_s2: accel,
        })
    }
}
