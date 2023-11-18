use rand::thread_rng;
use rand_distr::{Distribution, StandardNormal};
use rayon::prelude::*;

/* まだ足していく */

#[derive(Debug, Copy, Clone)]
pub struct CalcInput {
    pub zero_rate: f64,
    pub vol: f64,
    pub term_annu: f64,
    pub strike: f64,
    pub underlying: f64,
}

pub fn mc_bs_asian_call(input: &CalcInput, time_step: usize, num_path: usize) -> f64 {
    let CalcInput {
        zero_rate,
        vol,
        term_annu,
        strike,
        underlying,
    } = *input;
    let delta_t = term_annu / time_step as f64;

    let und_paths: Vec<Vec<f64>> = vec![vec![0.0; time_step]; num_path];
    let vals: Vec<f64> = und_paths
        .into_par_iter()
        .map(|mut und_path| -> f64 {
            und_path[0] = underlying;
            for i in 1..time_step {
                let norm_rand: f64 = StandardNormal.sample(&mut thread_rng());
                und_path[i] = und_path[i - 1]
                    * ((zero_rate - 0.5 * vol.powi(2)) * delta_t
                        + vol * delta_t.sqrt() * norm_rand)
                        .exp();
            }
            (-zero_rate * term_annu).exp()
                * (und_path.par_iter().sum::<f64>() / time_step as f64 - strike).max(0.0)
        })
        .collect();

    vals.par_iter().sum::<f64>() / num_path as f64
}
