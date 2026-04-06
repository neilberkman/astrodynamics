use nalgebra::Vector3;
use crate::state::CartesianState;
use crate::propagator::api::PropagationContext;
use crate::error::PropagationError;

pub trait ForceModel: Send + Sync {
    fn acceleration(
        &self,
        state: &CartesianState,
        ctx: &PropagationContext,
    ) -> Result<Vector3<f64>, PropagationError>;
}
