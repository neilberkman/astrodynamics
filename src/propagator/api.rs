#[derive(Debug, Clone, Default)]
pub struct PropagationContext {
    // For future expansion: frame, atmosphere model, etc.
}

pub struct IntegratorOptions {
    pub abs_tol: f64,
    pub rel_tol: f64,
    pub min_step: f64,
    pub max_step: f64,
    pub initial_step: f64,
    pub max_steps: u32,
    pub dense_output: bool,
}

impl Default for IntegratorOptions {
    fn default() -> Self {
        Self {
            abs_tol: 1e-9,
            rel_tol: 1e-12,
            min_step: 1e-6,
            max_step: 3600.0,
            initial_step: 60.0,
            max_steps: 1_000_000,
            dense_output: false,
        }
    }
}
