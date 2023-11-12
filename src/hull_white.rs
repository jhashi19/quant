mod distribution;

use distribution::moro_inverse_cdf_std_normal;

pub fn run() {
    for i in 0..1000 {
        let x = i as f64 / 1000.0;
        let val = moro_inverse_cdf_std_normal(x);
        println!("{}: {}", x, val);
    }
}
