use ndarray::{Array1, Array2, ArrayBase, Dim, OwnedRepr};
use ndarray_linalg::{FactorizeInto, Scalar, Solve};

pub fn levenberg_marquardt(
    f: &dyn Fn(&f64, &Vec<f64>) -> f64,
    args: &Vec<f64>,           // fitさせるときに固定されるパラメータ
    params: &mut Vec<f64>,     // firさせるときに変動させるパラメータの初期値
    target_vals: &Array1<f64>, // argsごとの関数の値の収束先
    max_iter: usize,           // イテレーション回数の最大
    threshold: f64,            // ループ終了のための閾値
) {
    let mut square_diff_new = 0.0; //残差二乗和(更新後)
    let mut square_diff_old = 0.0; //残差二乗和(更新前)
    let mut lamda = 0.001; //fudge factor
    let mut func_vals: Array1<f64> = Array1::zeros(target_vals.len()); // argsごとの関数の値
    let mut residuals: Array1<f64> = Array1::zeros(target_vals.len()); // 残差ベクトル
    let mut num_it = 0; //イテレーション回数

    while num_it < max_iter {
        for (idx, arg) in args.iter().enumerate() {
            func_vals[idx] = f(arg, &params);
        }
        residuals = &func_vals - target_vals;

        let jacobian = jacobian(f, &args, &params);
        // println!("jacobian: \n{:?}", jacobian);
        // println!(
        //     "jacobian.t().dot(&jacobian): \n{:?}",
        //     jacobian.t().dot(&jacobian)
        // );
        let lm_mat: ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>> =
            jacobian.t().dot(&jacobian) + lamda * Array2::eye(params.len());
        // println!("lm_mat: \n{:?}", lm_mat);
        let lm_mat = lm_mat.factorize_into().unwrap();
        let delta_params = lm_mat.solve_into(-jacobian.t().dot(&residuals)).unwrap();

        for i in 0..delta_params.len() {
            params[i] += delta_params[i];
        }

        square_diff_new = residuals.iter().map(|diff| diff.powi(2)).sum();

        if (square_diff_new - square_diff_old).abs() < threshold {
            break;
        }
        lamda *= if square_diff_new > square_diff_old {
            10.0
        } else {
            0.1
        };
        square_diff_old = square_diff_new;
        num_it += 1;
    }
    println!("params:{:?}", params);
    println!("num_it:{}", num_it);
}

fn jacobian(f: &dyn Fn(&f64, &Vec<f64>) -> f64, args: &Vec<f64>, params: &Vec<f64>) -> Array2<f64> {
    let mut jac: Array2<f64> = Array2::zeros((args.len(), params.len()));
    for i in 0..args.len() {
        for j in 0..params.len() {
            // TODO パラメータを動かしていない純粋な関数の値を何度も計算しているので、要改善
            jac[[i, j]] = numerical_difference(f, &args[i], params, j);
        }
    }
    jac
}

fn numerical_difference(
    f: &dyn Fn(&f64, &Vec<f64>) -> f64,
    arg: &f64,
    params: &Vec<f64>,
    idx: usize,
) -> f64 {
    let mut h = 0.0;
    if params[idx] == 0.0 {
        h = 0.0001;
    } else {
        h = params[idx].abs() * 0.001;
    }

    let mut incre_params = params.clone();
    incre_params[idx] += h;

    // let mut decre_params = params.clone();
    // decre_params[idx] -= h;

    // (f(arg, &incre_params) - f(arg, &decre_params)) / h
    (f(arg, &incre_params) - f(arg, params)) / h
}

// fn gradient(
//     f: &Box<dyn Fn(&f64, &Array1<f64>) -> f64>,
//     arg: &f64,
//     params: &Array1<f64>,
// ) -> Vec<f64> {
//     let mut grad = vec![0.0; params.len()];
//     for idx in 0..grad.len() {
//         grad[idx] = numerical_difference(f, arg, params, idx);
//     }
//     grad
// }
