pub mod rk4;
pub mod dp54;
pub mod tableau;

pub use rk4::RK4;
pub use dp54::DP54;

use crate::state::{CartesianState, StateDerivative};
use crate::propagator::api::{PropagationContext, IntegratorOptions};
use crate::propagator::result::PropagationResult;
use crate::error::PropagationError;

pub trait DynamicsModel: Send + Sync {
    fn derivative(
        &self,
        state: &CartesianState,
        ctx: &PropagationContext,
    ) -> Result<StateDerivative, PropagationError>;
}

pub trait Integrator: Send + Sync {
    fn propagate(
        &self,
        initial: CartesianState,
        t_end_seconds: f64,
        rhs: &dyn DynamicsModel,
        ctx: &PropagationContext,
        opts: &IntegratorOptions,
    ) -> Result<PropagationResult, PropagationError>;
}
