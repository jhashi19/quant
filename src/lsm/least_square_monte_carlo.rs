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

    // 2次元Vecで並列にパスを生成。いくつか試した中ではこれが一番速い。
    let und_paths: Vec<Vec<f64>> = vec![vec![0.0; time_step + 1]; NUM_PATH];
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

    // Array1を使ってVecを要素としてpar_map_inplaceで並列にパスを生成。
    // let mut und_paths: Array1<Vec<f64>> = Array1::from_elem(NUM_PATH, vec![0.0; time_step + 1]);
    // und_paths.par_map_inplace(|und_path: &mut Vec<f64>| {
    //     und_path[0] = underlying;
    //     for i in 0..time_step {
    //         let norm_rand: f64 = StandardNormal.sample(&mut thread_rng());
    //         und_path[i + 1] = und_path[i]
    //             * ((zero_rate - 0.5 * vol.powi(2)) * delta_t + vol * delta_t.sqrt() * norm_rand)
    //                 .exp();
    //     }
    // });

    let payoffs: Vec<f64> = vec![0.0; NUM_PATH];
    // for i in 0..NUM_PATH {
    //     payoffs[i] = (strike - und_paths[i][time_step]).max(0.0);
    // }
    let mut payoffs: Vec<f64> = payoffs
        .par_iter()
        .enumerate()
        .map(|(idx, _)| (strike - &und_paths[idx][time_step]).max(0.0))
        .collect();

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
        let continuation_vals = laguerre.dot(&reg_coeffs);
        for (itm_path_idx, path) in itm_paths.iter().enumerate() {
            // itm_pathsはITMのみのため、0との比較は不要
            let intrinsic_val = strike - und_paths[*path][step];
            if intrinsic_val > continuation_vals[itm_path_idx] {
                payoffs[*path] = intrinsic_val;
                not_excercise[*path] = false;
            }
        }

        for i in 0..NUM_PATH {
            if not_excercise[i] {
                payoffs[i] = df_one_step * payoffs[i];
            }
        }

        // 回帰をITMに絞らないパターン。ITMのパスのインデックスを抽出する必要がなくなるが、遅くなる。
        // let laguerre: Array2<f64> = Array2::from_shape_fn((NUM_PATH, 4), |(i, j)| {
        //     laguerre_polynomial(und_paths[i][step])[j]
        // });

        // let explained: Array1<f64> = Array1::from_shape_fn(NUM_PATH, |i| df_one_step * &payoffs[i]);
        // let reg_coeffs: Array1<f64> =
        //     ((((laguerre.t()).dot(&laguerre)).inv().unwrap()).dot(&laguerre.t())).dot(&explained);

        // let conti_vals = laguerre.dot(&reg_coeffs);

        // for i in 0..NUM_PATH {
        //     let intrinsic_val = strike - und_paths[i][step];
        //     if intrinsic_val > conti_vals[i] {
        //         payoffs[i] = intrinsic_val;
        //     } else {
        //         payoffs[i] = df_one_step * payoffs[i];
        //     }
        // }
    }
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
