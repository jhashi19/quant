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

    pub fn discount_bond_option(
        a: f64,
        sigma: f64,
        mat_u: f64,
        mat_o: f64,
        strike: f64,
        op_type: analysis::OptionType,
    ) -> f64 {
        let vol = analysis::dbo_vol(a, sigma, mat_u, mat_o);
        analysis::discount_bond_option(mat_u, mat_o, strike, vol, op_type)
    }

    pub fn capfloorlet(
        a: f64,
        sigma: f64,
        date_s: f64,
        date_e: f64,
        strike: f64,
        cf_type: analysis::CapFloorType,
        curve: curve::Curve,
    ) -> f64 {
        let vol = analysis::capfloorlet_vol(a, sigma, date_s, date_e);
        analysis::capfloorlet(date_s, date_e, strike, vol, cf_type, curve)
    }

    pub fn capfloor(
        dates: &Vec<f64>,
        vols: &Vec<f64>,
        strike: f64,
        cf_type: analysis::CapFloorType,
        curve: curve::Curve,
    ) -> f64 {
        analysis::capfloor(dates, vols, strike, cf_type, curve)
    }

    pub fn swaption(
        a: f64,
        sigma: f64,
        mat_op: f64,
        strike: f64,
        swap_dates: &Vec<f64>,
        op_type: analysis::SwaptionType,
        curve: curve::Curve,
    ) -> f64 {
        let vol = analysis::swaption_vol(a, sigma, mat_op, swap_dates, strike, curve);
        analysis::swaption(swap_dates, strike, vol, op_type, curve)
    }
}
