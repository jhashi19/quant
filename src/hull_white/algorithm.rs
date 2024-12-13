use ndarray::{Array1, Array2, ArrayBase, Dim, OwnedRepr};
use ndarray_linalg::{FactorizeInto, Solve};

pub struct Newton {
    func: fn(x: f64) -> f64, // ターゲットとなる1次元の関数。 =0の形であらわした左辺。
    func_deriv: fn(x: f64) -> f64, // funcの導関数
    x_max: f64,              // funcの表す方程式の解が取りうる値の範囲の上限
    x_min: f64,              // funcの表す方程式の解が取りうる値の範囲の下限
    threshold: f64,          // 収束判定の閾値
}

impl Newton {
    pub fn new(
        func: fn(x: f64) -> f64,
        func_deriv: fn(x: f64) -> f64,
        x_max: f64,
        x_min: f64,
        threshold: f64,
    ) -> Self {
        Newton {
            func,
            func_deriv,
            x_max,
            x_min,
            threshold,
        }
    }

    pub fn find_root(&self) -> f64 {
        let max_iter = 100; // 最大イテレーション回数
        let mut root = 0.5 * (self.x_max + self.x_min); // 根の初期値は根が存在する区間の中間
        for _ in 0..max_iter {
            let f = (self.func)(root);
            let df = (self.func_deriv)(root);
            let dx = f / df;
            root -= dx;
            if (self.x_max - root) * (root - self.x_min) < 0.0 {
                panic!("Root is out of range [{},{}]", self.x_min, self.x_max);
            }
            if dx.abs() < self.threshold {
                return root;
            }
        }
        panic!("Maximun number of iterations exceeded");
    }

    pub fn find_root_safe(&self) -> f64 {
        let max_iter = 100;
        let mut xf_high: f64; // x_maxとx_minのうち、funcが大きい値をとる方
        let mut xf_low: f64; // x_maxとx_minのうち、funcが小さい値をとる方
        let fr = (self.func)(self.x_max);
        let fl = (self.func)(self.x_min);
        if (fl > 0.0 && fr > 0.0) || (fl < 0.0 && fr < 0.0) {
            panic!("Root must be bracketed in [{},{}]", self.x_min, self.x_max);
        }
        if fl == 0.0 {
            return self.x_min;
        }
        if fr == 0.0 {
            return self.x_max;
        }
        if fl < 0.0 {
            xf_high = self.x_max;
            xf_low = self.x_min;
        } else {
            xf_high = self.x_min;
            xf_low = self.x_max;
        }
        let mut root = 0.5 * (self.x_max + self.x_min);
        let mut dxold = self.x_max - self.x_min;
        let mut dx = dxold;
        let mut f = (self.func)(root);
        let mut df = (self.func_deriv)(root);
        for _ in 0..max_iter {
            if ((root - xf_high) * df - f) * ((root - xf_low) * df - f) > 0.0
                || (2.0 * f).abs() > (dxold * df).abs()
            {
                dxold = dx;
                dx = 0.5 * (xf_high - xf_low);
                root = xf_low + dx;
            } else {
                dxold = dx;
                dx = f / df;
                root -= dx;
            }
            if dx.abs() < self.threshold {
                return root;
            }
            f = (self.func)(root);
            df = (self.func_deriv)(root);
            if f < 0.0 {
                xf_low = root;
            } else {
                xf_high = root;
            }
        }
        panic!("Maximun number of iterations exceeded");
    }
}

// struct では実装しづらいため関数で実装。
pub fn LevenbergMarquardt(
    func: &dyn Fn(f64, &Vec<f64>) -> f64, // キャリブレーションする対象の関数
    func_args: &Vec<f64>, // 関数の引数のうち、キャリブレーションに使用する独立変数となる変数以外（調整パラメータを含む）
    are_adjustings: Vec<bool>, // 関数の引数のうち、キャリブレーションで調整するパラメータのインデックスをtrueとするベクタ
    derivative_funcs: &Vec<&dyn Fn(f64, &Vec<f64>) -> f64>, // 各調整パラメータごとの、関数の偏導関数のベクタ
    independent_vars: &Vec<f64>, // キャリブレーションに使用する独立変数のベクタ
    dependent_vars: &Vec<f64>,   // キャリブレーションに使用する従属変数のベクタ
) -> Vec<f64> {
    let mut squared_sum = 0.0; //誤差の2乗和
    let mut new_squared_sum: f64; //更新後の誤差の2乗和
    let mut errors: Array1<f64>; // 従属変数と関数の値の差分
    let mut lambda = 0.001; //damping parameter
    let mut iteration_num = 0; // イテレーション回数（ログに出力するためにカウントする）
    let mut adjusted_args = Vec::with_capacity(func_args.len());
    let max_iter = 1000; // 最大イテレーション回数。これを超えると処理は異常終了とする。
    let threshold = 10e-7; // 誤差のイテレーションごとの変化が小さくなった場合にループを終了する閾値

    while iteration_num < max_iter {
        errors = calc_error(func, independent_vars, dependent_vars, func_args);

        let jacobian = jacobian(derivative_funcs, func_args, independent_vars);

        adjusted_args =
            update_adjusted_params(jacobian, &errors, func_args, &are_adjustings, lambda);

        new_squared_sum = errors.iter().map(|diff| diff.powi(2)).sum();

        if (new_squared_sum - squared_sum).abs() < threshold {
            println!("params:{:?}", adjusted_args);
            println!("iteration_num:{}", iteration_num);
            return adjusted_args;
        }

        lambda = update_lambda(squared_sum, new_squared_sum, lambda);
        squared_sum = new_squared_sum;
        iteration_num += 1;
    }
    println!("params:{:?}", adjusted_args);
    println!("iteration_num:{}", iteration_num);
    adjusted_args
}

/// 引数のベクタをArray1型に変換して返します。
fn vec_to_array(vec: &Vec<f64>) -> Array1<f64> {
    Array1::from_shape_fn(vec.len(), |i| vec[i])
}

/// 関数の値と従属変数の差異を返します。
fn calc_error(
    func: &dyn Fn(f64, &Vec<f64>) -> f64,
    independent_vars: &Vec<f64>,
    dependent_vars: &Vec<f64>,
    func_args: &Vec<f64>,
) -> Array1<f64> {
    // 独立変数ごとの関数の値（パラメータを調整するごとに更新）
    let mut func_vals: Array1<f64> = Array1::zeros(independent_vars.len());
    for i in 0..independent_vars.len() {
        func_vals[i] = func(independent_vars[i], func_args);
    }
    &func_vals - &vec_to_array(dependent_vars)
}

fn update_adjusted_params(
    jacobian: Array2<f64>,
    errors: &Array1<f64>,
    func_args: &Vec<f64>,
    are_adjustings: &Vec<bool>,
    lambda: f64,
) -> Vec<f64> {
    let mut updated_args = func_args.clone();

    let modified_jtj: ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>> =
        jacobian.t().dot(&jacobian) + lambda * Array2::eye(jacobian.shape()[1]);

    // let jtj = jacobian.t().dot(&jacobian);
    // let modified_jtj = jacobian.t().dot(&jacobian) + lamda * modified_jtj.diag();

    let modified_jtj = modified_jtj.factorize_into().unwrap();
    let delta_args = modified_jtj.solve_into(-jacobian.t().dot(errors)).unwrap();

    // TODO filterとかmapとかzipとか使って書けないか検討
    let mut delta_index = 0;
    for i in 0..are_adjustings.len() {
        if are_adjustings[i] {
            updated_args[i] += delta_args[delta_index];
            delta_index += 1;
        }
    }
    updated_args
}

fn update_lambda(squared_sum: f64, new_squared_sum: f64, lambda: f64) -> f64 {
    let lambda_up = 2.0; // damping parameter を大きくするときの掛ける値(delayed gratification)
    let lambda_down = 1.0 / 3.0; // damping parameter を小さくするときの掛ける値(delayed gratification)

    if new_squared_sum > squared_sum {
        lambda * lambda_up
    } else {
        lambda * lambda_down
    }
}

/// ヤコビアンを計算して返します。
fn jacobian(
    derivative_funcs: &Vec<&dyn Fn(f64, &Vec<f64>) -> f64>, // 各調整パラメータごとの偏導関数のベクタ
    args: &Vec<f64>,                                        // 偏導関数の引数
    independent_vars: &Vec<f64>,                            // 調整に使用する独立変数のベクタ
) -> Array2<f64> {
    let row_num = independent_vars.len();
    let column_num = derivative_funcs.len();
    let mut jacobian: Array2<f64> = Array2::zeros((row_num, column_num));
    for i in 0..row_num {
        for j in 0..derivative_funcs.len() {
            jacobian[[i, j]] = derivative_funcs[j](independent_vars[i], args);
        }
    }
    jacobian
}

/// 引数の関数を引数の1つで数値微分により偏微分した値を返します。
pub fn numerical_difference(f: &dyn Fn(&Vec<f64>) -> f64, args: &Vec<f64>, idx: usize) -> f64 {
    let h;
    if args[idx] == 0.0 {
        h = 0.0001;
    } else {
        h = args[idx].abs() * 0.001;
    }

    let mut incremented_args = args.clone();
    incremented_args[idx] += h;

    (f(&incremented_args) - f(args)) / h
}

#[cfg(test)]
mod tests {
    use super::*;

    fn func(x: f64) -> f64 {
        x.powi(2) - 2.0
    }
    fn func_deriv(x: f64) -> f64 {
        2.0 * x
    }

    #[test]
    fn test_newton_find_root() {
        let ins = Newton::new(func, func_deriv, 2.0, 0.0, 1.0e-7);
        let root = ins.find_root();
        assert!((root - 2_f64.powf(0.5)).abs() < 1.0e-10);
    }

    #[test]
    fn test_newton_find_root_safe() {
        let ins = Newton::new(func, func_deriv, 2.0, 0.0, 1.0e-7);
        let root = ins.find_root_safe();
        assert!((root - 2_f64.powf(0.5)).abs() < 1.0e-10);
    }
}
