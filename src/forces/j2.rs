use nalgebra::Vector3;
use crate::state::CartesianState;
use crate::propagator::api::PropagationContext;
use crate::error::PropagationError;
use crate::forces::r#trait::ForceModel;
use crate::constants::{MU_EARTH, RE_EARTH, J2_EARTH};

pub struct J2Gravity {
    pub mu: f64,
    pub re: f64,
    pub j2: f64,
}

impl Default for J2Gravity {
    fn default() -> Self {
        Self {
            mu: MU_EARTH,
            re: RE_EARTH,
            j2: J2_EARTH,
        }
    }
}

impl ForceModel for J2Gravity {
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
        
        let f = 1.5 * self.j2 * (self.mu / r_mag2) * (self.re / r_mag).powi(2);
        let z_r2 = (state.position_km.z * state.position_km.z) / r_mag2;
        
        Ok(Vector3::new(
            f * (state.position_km.x / r_mag) * (5.0 * z_r2 - 1.0),
            f * (state.position_km.y / r_mag) * (5.0 * z_r2 - 1.0),
            f * (state.position_km.z / r_mag) * (5.0 * z_r2 - 3.0),
        ))
    }
}
