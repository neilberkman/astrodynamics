use nalgebra::Vector3;
use crate::state::CartesianState;
use crate::propagator::api::PropagationContext;
use crate::error::PropagationError;
use crate::forces::r#trait::ForceModel;
use crate::constants::MU_EARTH;

pub struct TwoBodyGravity {
    pub mu: f64,
}

impl Default for TwoBodyGravity {
    fn default() -> Self {
        Self { mu: MU_EARTH }
    }
}

impl ForceModel for TwoBodyGravity {
    fn acceleration(
        &self,
        state: &CartesianState,
        _ctx: &PropagationContext,
    ) -> Result<Vector3<f64>, PropagationError> {
        let r_mag2 = state.position_km.norm_squared();
        if r_mag2 == 0.0 {
            return Err(PropagationError::NumericalFailure("Zero position magnitude".to_string()));
        }
        let r_mag = r_mag2.sqrt();
        Ok(state.position_km * (-self.mu / (r_mag2 * r_mag)))
    }
}
