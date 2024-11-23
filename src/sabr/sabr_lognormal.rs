// Hagan et al. 2002
pub fn sabr_lognormal(
    strike: f64,
    fwd: f64,
    term: f64,
    beta: f64,
    alpha: f64,
    rho: f64,
    nu: f64,
) -> f64 {
    let one_beta = 1.0 - beta;
    if (fwd - strike).abs() > 1e-5 {
        // ATMでない
        let fwd_mul_strike = fwd * strike;
        let fwd_div_strike = fwd / strike;
        let z = nu / alpha * fwd_mul_strike.powf(one_beta * 0.5) * (fwd_div_strike).ln();
        let x_z = (((1.0 - 2.0 * rho * z + z.powi(2)).sqrt() + z - rho) / (1.0 - rho)).ln();
        alpha
            * (1.0
                + term
                    * ((one_beta * alpha).powi(2) / (24.0 * fwd_mul_strike.powf(one_beta))
                        + rho * beta * nu * alpha / (4.0 * fwd_mul_strike.powf(one_beta * 0.5))
                        + (2.0 - 3.0 * rho.powi(2)) * nu.powi(2) / 24.0))
            / (fwd_mul_strike.powf(one_beta * 0.5)
                * (1.0
                    + (one_beta * fwd_div_strike.ln()).powi(2) / 24.0
                    + (one_beta * fwd_div_strike.ln()).powi(4) / 1920.0))
            * z
            / x_z
    } else {
        //ATM
        alpha / fwd.powf(one_beta)
            * (1.0
                + term
                    * (((one_beta * alpha).powi(2) / (24.0 * fwd.powf(2.0 * one_beta)))
                        + alpha * beta * nu * rho / (4.0 * fwd.powf(one_beta))
                        + (2.0 - 3.0 * rho.powi(2)) * nu.powi(2) / 24.0))
    }
}
