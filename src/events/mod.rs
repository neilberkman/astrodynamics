pub mod r#trait;

#[derive(Debug, Clone)]
pub struct DetectedEvent {
    pub epoch_tdb_seconds: f64,
    pub name: String,
    // Additional fields as needed
}
