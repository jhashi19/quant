mod algorithm;
mod analysis;
mod data;
mod distribution;
mod hw_lib;
mod tree;

use distribution::inverse_std_normal_cdf;

pub fn run() {
    for i in 0..1000 {
        let x = i as f64 / 1000.0;
        let val = inverse_std_normal_cdf(x);
        println!("{}: {}", x, val);
    }
}
