use super::analysis;
use super::data;
use super::tree;

pub struct HullWhite {
    pub a: f64,
    pub sigma: f64,
}

impl HullWhite {
    pub fn new(a: f64, sigma: f64) -> Self {
        HullWhite { a, sigma }
    }

    pub fn discount_bond_option(
        &self,
        mat_u: &f64,
        mat_o: &f64,
        strike: &f64,
        op_type: analysis::OptionType,
    ) -> f64 {
        let vol = analysis::dbo_sigma(&self.a, &self.sigma, mat_u, mat_o);
        analysis::discount_bond_option(mat_u, mat_o, strike, &vol, op_type)
    }

    pub fn capfloorlet(
        &self,
        date_s: &f64,
        date_e: &f64,
        strike: &f64,
        cf_type: analysis::CapFloorType,
        curve: data::Curve,
    ) -> f64 {
        let vol = analysis::capfloorlet_sigma(&self.a, &self.sigma, date_s, date_e);
        analysis::capfloorlet(date_s, date_e, strike, &vol, cf_type, curve)
    }

    pub fn capfloor(
        &self,
        dates: &Vec<f64>,
        vols: &Vec<f64>,
        strike: &f64,
        cf_type: analysis::CapFloorType,
        curve: data::Curve,
    ) -> f64 {
        analysis::capfloor(dates, vols, strike, cf_type, curve)
    }

    pub fn swaption(&self) -> f64 {
        // analysis::swaption()
        0.0
    }
}
