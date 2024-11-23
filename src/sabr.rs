mod lm;
mod sabr_lognormal;

use lm::levenberg_marquardt;
use ndarray::arr1;
use sabr_lognormal::sabr_lognormal;

pub fn run() {
    let fwd = 0.03571;
    let term = 10.0;
    let beta = 0.5;
    let strikes = vec![
        fwd - 0.02,
        fwd - 0.01,
        fwd - 0.005,
        fwd,
        fwd + 0.005,
        fwd + 0.01,
        fwd + 0.02,
    ];
    let mut params = vec![0.01, 0.01, 0.5];
    let target_vals = arr1(&[0.3215, 0.248, 0.2222, 0.2040, 0.1923, 0.1867, 0.1887]);
    let max_iter: usize = 1000;
    let threshold = 1e-7;
    let rap_fn = |arg: &f64, params: &Vec<f64>| -> f64 {
        sabr_lognormal(*arg, fwd, term, beta, params[0], params[1], params[2])
    };
    levenberg_marquardt(
        &rap_fn,
        &strikes,
        &mut params,
        &target_vals,
        max_iter,
        threshold,
    );
    for i in 0..strikes.len() {
        println!(
            "(sabr, target): {}, {}",
            rap_fn(&strikes[i], &params),
            target_vals[i]
        );
    }
}
