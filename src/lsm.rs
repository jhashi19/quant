mod least_square_monte_carlo;

use least_square_monte_carlo::{longstaff_schwartz_american_put, CalcInput};
use std::time::Instant;

pub fn run() {
    let input = CalcInput {
        underlying: 100.0,
        strike: 100.0,
        vol: 0.2,
        zero_rate: 0.05,
        term_annu: 1.0,
    };
    let start = Instant::now();
    let opt_price = longstaff_schwartz_american_put(&input, 100, 10000);
    let end = start.elapsed();
    println!("(lsm) time:{}s", end.as_secs_f64());
    println!("(lsm) american option price: {}", opt_price);
}
