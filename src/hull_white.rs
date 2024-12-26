mod analysis;
mod curve;
mod data;
mod hw_lib;
mod interpolation;
mod math;
mod node;
mod optimization;
mod tree;

use math::inverse_std_normal_cdf;

pub fn run() {
    for i in 0..1000 {
        let x = i as f64 / 1000.0;
        let val = inverse_std_normal_cdf(x);
        println!("{}: {}", x, val);
    }
}
