#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Time {
    pub seconds_since_j2000: f64,
}

impl Time {
    pub fn new(seconds_since_j2000: f64) -> Self {
        Self { seconds_since_j2000 }
    }

    pub fn tdb(&self) -> f64 {
        self.seconds_since_j2000
    }
}
