use crate::state::CartesianState;
use crate::events::DetectedEvent;
use crate::propagator::dense_output::DenseOutput;

#[derive(Debug, Clone)]
pub struct PropagationPoint {
    pub epoch_tdb_seconds: f64,
    pub position_km: [f64; 3],
    pub velocity_km_s: [f64; 3],
}

#[derive(Debug, Clone, Default)]
pub struct PropagationStats {
    pub accepted_steps: u32,
    pub rejected_steps: u32,
    pub evaluations: u32,
}

#[derive(Debug, Clone)]
pub struct PropagationResult {
    pub final_state: CartesianState,
    pub points: Vec<PropagationPoint>,
    pub events: Vec<DetectedEvent>,
    pub stats: PropagationStats,
    pub dense: Option<DenseOutput>,
}
