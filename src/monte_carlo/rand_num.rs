use rand::prelude::*;
use rand_distr::StandardNormal;

pub fn normal_rand() {
    let num = 1000000;
    let rng = thread_rng();
    let norm_rands: Vec<f64> = StandardNormal.sample_iter(rng).take(num).collect();
    let mean = norm_rands.iter().sum::<f64>() / num as f64;
    let vol = (norm_rands.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / num as f64).sqrt();
    println!("mean: {}", mean);
    println!("vol: {}", vol);
}
