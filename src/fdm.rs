use finite_difference_method::CalcInput;
use std::time::Instant;

mod finite_difference_method;

pub fn run() {
    let input = CalcInput {
        underlying: 62.0,
        strike: 60.0,
        vol: 0.2,
        zero_rate: 0.1,
        term_annu: 5.0 / 12.0,
    };
    let start = Instant::now();
    let value = finite_difference_method::explicit_fdm_bs(&input, 120.0, 2400, 96000);
    let end = start.elapsed();
    println!("(explicit_fdm_bs) time:{}s", end.as_secs_f64());
    println!("(explicit_fdm_bs) european option price: {}", value);
}
