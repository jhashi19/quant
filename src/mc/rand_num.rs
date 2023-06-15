use ndarray::Array1;
use rand::prelude::*;
use rand_distr::StandardNormal;
use rayon::prelude::*;

pub fn normal_rand() {
    const NUM: usize = 1000000;
    let mut norm_rands: Array1<f64> = Array1::zeros(NUM);
    norm_rands.par_mapv_inplace(|_| StandardNormal.sample(&mut thread_rng()));
    let mean = norm_rands.par_iter().sum::<f64>() / NUM as f64;
    let vol = (norm_rands
        .par_iter()
        .map(|x| (x - mean).powi(2))
        .sum::<f64>()
        / NUM as f64)
        .sqrt();
    println!("mean: {}", mean);
    println!("vol: {}", vol);
}
