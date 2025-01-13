use super::analysis;
use super::curve;
use super::tree;

#[derive(Clone, Debug)]
pub struct HullWhite {
    pub a: Vec<f64>,
    pub a_interval: Vec<f64>,
    pub sigma: Vec<f64>,
    pub sigma_interval: Vec<f64>,
}

impl HullWhite {
    pub fn new(
        a: Vec<f64>,
        a_interval: Vec<f64>,
        sigma: Vec<f64>,
        sigma_interval: Vec<f64>,
    ) -> Self {
        HullWhite {
            a,
            a_interval,
            sigma,
            sigma_interval,
        }
    }
}
