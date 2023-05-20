use std::f64::consts::PI;

#[derive(Debug, Copy, Clone)]
pub struct CalcInput {
    pub zero_rate: f64,
    pub vol: f64,
    pub term_annu: f64,
    pub strike: f64,
    pub underlying: f64,
}

pub enum OptionType {
    Call,
    Put,
}
pub fn black_scholes(input: &CalcInput, option_type: OptionType) -> f64 {
    let CalcInput {
        zero_rate,
        vol,
        term_annu,
        strike,
        underlying,
    } = input;
    let d1 = ((underlying / strike).ln() + (zero_rate + 0.5 * vol.powi(2)) * term_annu)
        / (vol * term_annu.sqrt());
    let d2 = d1 - vol * term_annu.sqrt();
    let sgn = match option_type {
        OptionType::Call => 1.0,
        OptionType::Put => -1.0,
    };
    sgn * (underlying * norm_cdf_matic2016(sgn * d1)
        - strike * (-zero_rate * term_annu).exp() * norm_cdf_matic2016(sgn * d2))
}

// Matic et al.(2016)
pub fn norm_cdf_matic2016(x: f64) -> f64 {
    let gamma2 = -1.0 / 3.0 + 1.0 / PI;
    let gamma4 = 7.0 / 90.0 - 2.0 / (3.0 * PI) + 4.0 / (3.0 * PI.powi(2));
    let gamma6 = -1.0 / 70.0 + 4.0 / (15.0 * PI) - 4.0 / (3.0 * PI.powi(2)) + 2.0 / PI.powi(3);
    let gamma8 = 83.0 / 37800.0 - 76.0 / (945.0 * PI) + 34.0 / (45.0 * PI.powi(2))
        - 8.0 / (3.0 * PI.powi(3))
        + 16.0 / (5.0 * PI.powi(4));
    let gamma10 = -73.0 / 249480.0 + 283.0 / (14175.0 * PI) - 178.0 / (567.0 * PI.powi(2))
        + 88.0 / (45.0 * PI.powi(3))
        - 16.0 / (3.0 * PI.powi(4))
        + 16.0 / (3.0 * PI.powi(5));

    0.5 + x.signum()
        * 0.5
        * (1.0
            - (-2.0 * x.powi(2) / PI
                * (1.0
                    + gamma2 * x.powi(2)
                    + gamma4 * x.powi(4)
                    + gamma6 * x.powi(6)
                    + gamma8 * x.powi(8)
                    + gamma10 * x.powi(10)))
            .exp())
        .sqrt()
}
