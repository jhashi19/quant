mod monte_carlo;
use monte_carlo::{mc_bs_asian_call, CalcInput};
use std::time::Instant;

// mod rand_num;

pub fn run() {
    // rand_num::normal_rand();
    let input = CalcInput {
        underlying: 100.0,
        strike: 100.0,
        vol: 0.2,
        zero_rate: 0.05,
        term_annu: 1.0,
    };
    let start = Instant::now();
    let opt_price = mc_bs_asian_call(&input, 250, 10000);
    let end = start.elapsed();
    println!("(monte_carlo) time:{}s", end.as_secs_f64());
    println!("(monte_carlo) mc_bs_asian_call: {}", opt_price);
}
