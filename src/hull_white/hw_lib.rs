use super::analysis;
use super::tree;

pub struct HullWhite {
    pub a: f64,
    pub sigma: f64,
    theta: f64,
}

impl HullWhite {
    pub fn df(&self, t: &f64, matu: &f64) -> f64 {
        analysis::df(t, matu)
    }

    pub fn discount_bond_option(
        &self,
        t: &f64,
        mat_u: &f64,
        mat_o: &f64,
        strike: &f64,
        op_type: analysis::OptionType,
    ) -> f64 {
        analysis::discount_bond_option(&self.a, &self.sigma, t, mat_u, mat_o, strike, op_type)
    }

    pub fn capfloorlet(&self) -> f64 {
        analysis::capfloorlet()
    }

    pub fn capfloor(&self) -> f64 {
        analysis::capfloor()
    }

    pub fn swaption(&self) -> f64 {
        analysis::swaption()
    }
}
