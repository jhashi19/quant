use ndarray::{Array1, Array2};
use ndarray_linalg::Inverse;
use rand::thread_rng;
use rand_distr::{Distribution, StandardNormal};
use rayon::prelude::*;

#[derive(Debug, Copy, Clone)]
pub struct CalcInput {
    pub zero_rate: f64,
    pub vol: f64,
    pub term_annu: f64,
    pub strike: f64,
    pub underlying: f64,
}

pub fn longstaff_schwartz_american_put(input: &CalcInput, time_step: usize) -> f64 {
    let CalcInput {
        zero_rate,
        vol,
        term_annu,
        strike,
        underlying,
    } = *input;
    const NUM_PATH: usize = 10000;
    println!("longstaff_schwartz_american_put");
    println!("number of paths:{}", NUM_PATH);

    let delta_t = term_annu / time_step as f64;
    let df_one_step = (-zero_rate * delta_t).exp();

    // let mut und_paths: Array2<f64> = Array2::zeros((NUM_PATH, time_step + 1));
    let und_paths: Vec<Vec<f64>> = vec![vec![0.0; time_step + 1]; NUM_PATH];
    // let mut payoffs: Array1<f64> = Array1::zeros(NUM_PATH);
    // let payoffs_tmp: Vec<f64> = vec![0.0; NUM_PATH];
    let mut payoffs: Vec<f64> = vec![0.0; NUM_PATH];
    // let mut norm_rand: Array2<f64> = Array2::zeros((NUM_PATH, time_step));
    // norm_rand.par_mapv_inplace(|_| -> f64 { StandardNormal.sample(&mut thread_rng()) });
    let und_paths: Vec<Vec<f64>> = und_paths
        .into_par_iter()
        .map(|mut und_path| -> Vec<f64> {
            und_path[0] = underlying;
            for i in 0..time_step {
                let norm_rand: f64 = StandardNormal.sample(&mut thread_rng());
                und_path[i + 1] = und_path[i]
                    * ((zero_rate - 0.5 * vol.powi(2)) * delta_t
                        + vol * delta_t.sqrt() * norm_rand)
                        .exp();
            }
            und_path
        })
        .collect();
    for i in 0..NUM_PATH {
        payoffs[i] = (strike - und_paths[i][time_step]).max(0.0);
    }
    // let payoffs_tmp: Vec<f64> = payoffs_tmp
    //     .par_iter()
    //     .enumerate()
    //     .map(|(idx, _)| (strike - &und_paths[idx][time_step]).max(0.0))
    //     .collect();
    // let mut payoffs = &mut payoffs_tmp.clone();
    // drop(payoffs_tmp);

    // for i in 0..NUM_PATH {
    //     und_paths[[i, 0]] = underlying;
    //     for j in 1..=time_step {
    //         und_paths[[i, j]] = und_paths[[i, j - 1]]
    //             * ((zero_rate - 0.5 * vol.powi(2)) * delta_t
    //                 + vol * delta_t.sqrt() * norm_rand[[i, j - 1]])
    //             .exp();
    //     }
    //     payoffs[i] = (strike - und_paths[[i, time_step]]).max(0.0);
    // }

    for step in (1..time_step).rev() {
        let mut itm_paths: Vec<usize> = Vec::new();
        for path in 0..NUM_PATH {
            if strike - und_paths[path][step] > 0.0 {
                itm_paths.push(path);
            }
        }
        // 基底関数はラゲール多項式の0～3項を使用する。
        let laguerre: Array2<f64> = Array2::from_shape_fn((itm_paths.len(), 4), |(i, j)| {
            laguerre_polynomial(und_paths[itm_paths[i]][step])[j]
        });

        let explained: Array1<f64> =
            Array1::from_shape_fn(itm_paths.len(), |i| df_one_step * &payoffs[itm_paths[i]]);
        let reg_coeffs: Array1<f64> =
            ((((laguerre.t()).dot(&laguerre)).inv().unwrap()).dot(&laguerre.t())).dot(&explained);

        let mut not_excercise = [true; NUM_PATH];
        for (itm_path_idx, path) in itm_paths.iter().enumerate() {
            // itm_pathsはITMのみのため、0との比較は不要
            let intrinsic_val = strike - und_paths[*path][step];
            let mut continuation_val = 0.0;
            for i in 0..=3 {
                continuation_val += reg_coeffs[i] * laguerre[[itm_path_idx, i]];
            }
            if intrinsic_val > continuation_val {
                payoffs[*path] = intrinsic_val;
                not_excercise[*path] = false;
            }
        }

        for i in 0..NUM_PATH {
            if not_excercise[i] {
                payoffs[i] = df_one_step * payoffs[i];
            }
        }
    }
    // println!("{:?}", &payoffs);
    payoffs.par_iter().sum::<f64>() / NUM_PATH as f64 * df_one_step
}

// ラゲール多項式の0〜3までの項を返す
fn laguerre_polynomial(x: f64) -> [f64; 4] {
    let lag0 = 1.0;
    let lag1 = 1.0 - x;
    let lag2 = 1.0 - 2.0 * x + 0.5 * x.powi(2);
    let lag3 = 1.0 - 3.0 * x + 1.5 * x.powi(2) - x.powi(3) / 6.0;
    [lag0, lag1, lag2, lag3]
}
