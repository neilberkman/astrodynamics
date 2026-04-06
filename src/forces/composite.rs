use nalgebra::Vector3;
use crate::state::CartesianState;
use crate::propagator::api::PropagationContext;
use crate::error::PropagationError;
use crate::forces::r#trait::ForceModel;

pub struct CompositeForceModel {
    pub models: Vec<Box<dyn ForceModel>>,
}

impl CompositeForceModel {
    pub fn new() -> Self {
        Self { models: Vec::new() }
    }

    pub fn add(&mut self, model: Box<dyn ForceModel>) {
        self.models.push(model);
    }
}

impl ForceModel for CompositeForceModel {
    fn acceleration(
        &self,
        state: &CartesianState,
        ctx: &PropagationContext,
    ) -> Result<Vector3<f64>, PropagationError> {
        let mut accel = Vector3::zeros();
        for model in &self.models {
            accel += model.acceleration(state, ctx)?;
        }
        Ok(accel)
    }
}
