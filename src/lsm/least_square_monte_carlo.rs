use ndarray::{Array1, Array2};
use ndarray_linalg::Inverse;
use rand::prelude::ThreadRng;
use rand_distr::{DistIter, Distribution, StandardNormal};

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
    const ITER_NUM: usize = 10000;
    println!("longstaff_schwartz_ap");
    println!("iterate:{}", ITER_NUM);

    let delta_t = term_annu / time_step as f64;
    let df_one_step = (-zero_rate * delta_t).exp();

    let mut und_paths: Array2<f64> = Array2::from_elem((ITER_NUM, time_step + 1), 0.0);
    let mut payoffs: Array1<f64> = Array1::from_elem(ITER_NUM, 0.0);
    let mut rng = rand::thread_rng();
    for i in 0..ITER_NUM {
        let mut norm_rand: DistIter<StandardNormal, &mut ThreadRng, f64> =
            StandardNormal.sample_iter(&mut rng);
        und_paths[[i, 0]] = underlying;
        for j in 1..=time_step {
            und_paths[[i, j]] = und_paths[[i, j - 1]]
                * ((zero_rate - 0.5 * vol.powi(2)) * delta_t
                    + vol * delta_t.sqrt() * norm_rand.next().unwrap())
                .exp();
        }
        payoffs[i] = (strike - und_paths[[i, time_step]]).max(0.0);
    }

    for step in (1..time_step).rev() {
        let mut itm_paths: Vec<usize> = Vec::new();
        for path in 0..ITER_NUM {
            if strike - und_paths[[path, step]] > 0.0 {
                itm_paths.push(path);
            }
        }
        // 基底関数はラゲール多項式の0～3項を使用する。
        let laguerre: Array2<f64> = Array2::from_shape_fn((itm_paths.len(), 4), |(i, j)| {
            laguerre_polynomial(und_paths[[itm_paths[i], step]])[j]
        });

        let explained: Array1<f64> =
            Array1::from_shape_fn(itm_paths.len(), |i| df_one_step * payoffs[itm_paths[i]]);
        let reg_coeffs =
            ((((laguerre.t()).dot(&laguerre)).inv().unwrap()).dot(&laguerre.t())).dot(&explained);

        let mut not_excercise = [true; ITER_NUM];
        for (itm_path_idx, path) in itm_paths.iter().enumerate() {
            // itm_pathsはITMのみのため、0との比較は不要
            let intrinsic_val = strike - und_paths[[*path, step]];
            let mut continuation_val = 0.0;
            for i in 0..=3 {
                continuation_val += reg_coeffs[i] * laguerre[[itm_path_idx, i]];
            }
            if intrinsic_val > continuation_val {
                payoffs[*path] = intrinsic_val;
                not_excercise[*path] = false;
            }
        }
        for i in 0..ITER_NUM {
            if not_excercise[i] {
                payoffs[i] = df_one_step * payoffs[i];
            }
        }
    }
    payoffs.iter().sum::<f64>() / ITER_NUM as f64 * df_one_step
}

// ラゲール多項式の0〜3までの項を返す
fn laguerre_polynomial(x: f64) -> [f64; 4] {
    let lag0 = 1.0;
    let lag1 = 1.0 - x;
    let lag2 = 1.0 - 2.0 * x + 0.5 * x.powi(2);
    let lag3 = 1.0 - 3.0 * x + 1.5 * x.powi(2) - x.powi(3) / 6.0;
    [lag0, lag1, lag2, lag3]
}
