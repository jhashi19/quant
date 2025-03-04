use ndarray::{Array1, Array2};
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

// levenberg_marquardt で使用する関数の型
type FnTypeLM = fn(f64, &[f64]) -> f64;

// struct では実装しづらいため関数で実装。
pub fn levenberg_marquardt<F, G>(
    func: F,                    // キャリブレーションする対象の関数
    func_args: Vec<f64>, // 関数の引数のうち、キャリブレーションに使用する独立変数となる変数以外（調整パラメータを含む）
    are_adjustings: Vec<bool>, // 関数の引数のうち、キャリブレーションで調整するパラメータのインデックスをtrueとするベクタ
    derivative_funcs: Vec<G>, // 各調整パラメータごとの、関数の偏導関数のベクタ（調整する対象のパラメータの偏微分のみ）
    independent_vars: Vec<f64>, // キャリブレーションに使用する独立変数のベクタ
    dependent_vars: Vec<f64>, // キャリブレーションに使用する従属変数のベクタ
) -> Vec<f64>
where
    F: Fn(f64, &[f64]) -> f64,
    G: Fn(f64, &[f64]) -> f64,
{
    let mut errors = calc_error(&func, &independent_vars, &dependent_vars, &func_args); // 従属変数と関数の値の差分
    let mut squared_sum: f64 = errors.iter().map(|diff| diff.powi(2)).sum(); //誤差の2乗和
    let mut new_squared_sum: f64; //更新後の誤差の2乗和
    let mut jacobian = calc_jacobian(&derivative_funcs, &func_args, &independent_vars); // ヤコビアン
    let mut params = func_args; // キャリブレーション対象パラメータ
    let mut new_params: Vec<f64>; // 更新されたキャリブレーション対象パラメータ
    let mut lambda = 0.001; //damping parameter
    const LAMBDA_UP: f64 = 2.0; // damping parameter を大きくするときの掛ける値(delayed gratification)
    const LAMBDA_DOWN: f64 = 1.0 / 3.0; // damping parameter を小さくするときの掛ける値(delayed gratification)
    const MAX_ITER: usize = 1000; // 最大イテレーション回数。これを超えると処理は異常終了とする。
    const THRESHOLD: f64 = 1e-7; // 誤差のイテレーションごとの変化が小さくなった場合にループを終了する閾値
    let mut iteration_count = 0; // イテレーション回数

    while iteration_count < MAX_ITER {
        new_params = update_params(&jacobian, &errors, &params, &are_adjustings, lambda);
        errors = calc_error(&func, &independent_vars, &dependent_vars, &new_params);
        new_squared_sum = errors.iter().map(|diff| diff.powi(2)).sum();
        if (new_squared_sum - squared_sum).abs() < THRESHOLD {
            println!("params:{:?}", new_params);
            println!("iteration_count:{}", iteration_count);
            return new_params;
        }

        if new_squared_sum < squared_sum {
            params = new_params;
            jacobian = calc_jacobian(&derivative_funcs, &params, &independent_vars);
            lambda *= LAMBDA_DOWN;
            squared_sum = new_squared_sum;
        } else {
            lambda *= LAMBDA_UP;
        }
        iteration_count += 1;
    }
    println!("params:{:?}", params);
    println!("iteration_count:{}", iteration_count);
    panic!("Maximun number of iterations exceeded");
}

/// 引数のベクタをArray1型に変換して返します。
fn vec_to_array(vec: &Vec<f64>) -> Array1<f64> {
    Array1::from_shape_fn(vec.len(), |i| vec[i])
}

/// 現時点のパラメータをもとに計算した関数の値と従属変数の値の差異を返します。
fn calc_error<F>(
    func: &F,
    independent_vars: &Vec<f64>,
    dependent_vars: &Vec<f64>,
    func_args: &Vec<f64>,
) -> Array1<f64>
where
    F: Fn(f64, &[f64]) -> f64,
{
    // 独立変数ごとの関数の値（パラメータを調整するごとに更新）
    let mut func_vals: Array1<f64> = Array1::zeros(independent_vars.len());
    for i in 0..independent_vars.len() {
        func_vals[i] = func(independent_vars[i], func_args);
    }
    &func_vals - &vec_to_array(dependent_vars)
}

/// パラメータを更新して返します。
fn update_params(
    jacobian: &Array2<f64>,
    errors: &Array1<f64>,
    func_args: &Vec<f64>,
    are_adjustings: &Vec<bool>,
    lambda: f64,
) -> Vec<f64> {
    let mut updated_args = func_args.clone();

    let modified_jtj: Array2<f64> =
        jacobian.t().dot(jacobian) + lambda * Array2::eye(jacobian.shape()[1]);

    // 以下は各対角成分に一律でlambdaを足すのではなく、jtjの対角成分のlambda倍を足すパターン
    // let jtj = jacobian.t().dot(jacobian);
    // let modified_jtj = &jtj + &(&Array2::from_diag(&jtj.diag()) * lambda);

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

/// ヤコビアンを計算して返します。
fn calc_jacobian<F>(
    derivative_funcs: &Vec<F>,   // 各調整パラメータごとの偏導関数のベクタ
    args: &Vec<f64>,             // 偏導関数の引数
    independent_vars: &Vec<f64>, // 調整に使用する独立変数のベクタ
) -> Array2<f64>
where
    F: Fn(f64, &[f64]) -> f64,
{
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
pub fn numerical_difference<F>(f: F, args: &Vec<f64>, idx: usize) -> f64
where
    F: Fn(&[f64]) -> f64,
{
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

/// levenberg_marquardt用に数値微分による偏微分を計算するクロージャを返します。
/// 返すクロージャの最初の引数は独立変数、2番目は調整パラメータのベクタです。
fn numerical_difference_for_lm<F>(f: F, index: usize) -> impl Fn(f64, &[f64]) -> f64
where
    F: Fn(f64, &[f64]) -> f64,
{
    move |x: f64, v: &[f64]| {
        let mut params = v.to_vec();
        params.insert(0, x);
        numerical_difference(|args| f(args[0], &args[1..].to_vec()), &params, index)
    }
}

/// levenberg_marquardtのderivative_funcs引数にそのまま設定できる形式で数値微分による偏微分関数を返します。
pub fn derivative_funcs_numerical_difference_for_lm<'a, F>(
    f: &'a F,
    index_num: usize,
) -> Vec<Box<dyn Fn(f64, &[f64]) -> f64 + 'a>>
where
    F: Fn(f64, &[f64]) -> f64 + Copy + 'a,
{
    let mut derivative_funcs: Vec<Box<dyn Fn(f64, &[f64]) -> f64 + 'a>> =
        Vec::with_capacity(index_num);
    for i in 1..=index_num {
        derivative_funcs.push(Box::new(numerical_difference_for_lm(f, i)));
    }
    derivative_funcs
}

#[cfg(test)]
mod tests {
    use super::*;

    fn func_test_newton(x: f64) -> f64 {
        x.powi(2) - 2.0
    }
    fn func_test_newton_deriv(x: f64) -> f64 {
        2.0 * x
    }

    #[test]
    fn test_newton_find_root() {
        let ins = Newton::new(func_test_newton, func_test_newton_deriv, 2.0, 0.0, 1.0e-7);
        let root = ins.find_root();
        assert!((root - 2_f64.powf(0.5)).abs() < 1.0e-10);
    }

    #[test]
    fn test_newton_find_root_safe() {
        let ins = Newton::new(func_test_newton, func_test_newton_deriv, 2.0, 0.0, 1.0e-7);
        let root = ins.find_root_safe();
        assert!((root - 2_f64.powf(0.5)).abs() < 1.0e-10);
    }

    #[test]
    fn test_levenberg_marquardt() {
        fn f(w: f64, v: &[f64]) -> f64 {
            let x = v[0];
            let y = v[1];
            let z = v[2];
            x + w * y + w.powi(2) * z
        }
        let f_deriv_x = |w: f64, v: &[f64]| -> f64 { 1.0 };
        // let f_deriv_y = |w: f64, v: &[f64]| -> f64 { w };
        let f_deriv_z = |w: f64, v: &[f64]| -> f64 { w.powi(2) };

        let func_args = vec![1.0, 1.9, 1.0];
        let are_adjustings = vec![true, false, true];
        let derivative_funcs: Vec<FnTypeLM> = vec![f_deriv_x, f_deriv_z];
        let independent_vars = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let dependent_vars = vec![-0.9, 1.9, 7.3, 13.8, 23.5];
        let actual_vec = levenberg_marquardt(
            &f,
            func_args,
            are_adjustings,
            derivative_funcs,
            independent_vars,
            dependent_vars,
        );
        let tolerance = 1e-10;
        assert!((actual_vec[0] - (-0.945517241379294)).abs() < tolerance);
        assert_eq!(actual_vec[1], 1.9);
        assert!((actual_vec[2] - 1.04425287356321).abs() < tolerance);
    }

    #[test]
    fn test_numerical_difference() {
        fn f(v: &[f64]) -> f64 {
            let x = v[0];
            let y = v[1];
            x.powi(2) + y.powi(3) - 2.0 * x * y
        }
        fn f_deriv_x(v: &[f64]) -> f64 {
            let x = v[0];
            let y = v[1];
            2.0 * x - 2.0 * y
        }
        fn f_deriv_y(v: &[f64]) -> f64 {
            let x = v[0];
            let y = v[1];
            3.0 * y.powi(2) - 2.0 * x
        }

        let args = vec![-0.24, 0.51];
        // 期待値との差の許容範囲はとりあえずで設定。
        // 引数の水準から決めようとしたが、対象の関数によるため。
        let tolerance = 1e-1;
        // let tolerance = args
        //     .iter()
        //     .map(|x: &f64| x.abs())
        //     .fold(0.0, |m, v| v.max(m))
        //     / 10.0;
        let actual_deriv_x = numerical_difference(f, &args, 0);
        let expected_deriv_x = f_deriv_x(&args);
        println!("actual_deriv_x: {}", actual_deriv_x);
        println!("expected_deriv_x: {}", expected_deriv_x);

        let actual_deriv_y = numerical_difference(f, &args, 1);
        let expected_deriv_y = f_deriv_y(&args);
        println!("actual_deriv_y: {}", actual_deriv_y);
        println!("expected_deriv_y: {}", expected_deriv_y);

        println!("tolerance: {}", tolerance);
        assert!((actual_deriv_x - expected_deriv_x).abs() < tolerance);
        assert!((actual_deriv_y - expected_deriv_y).abs() < tolerance);
    }

    #[test]
    fn test_levenberg_marquardt_numerical_difference() {
        fn f(w: f64, v: &[f64]) -> f64 {
            let x = v[0];
            let y = v[1];
            let z = v[2];
            x + w * y + w.powi(2) * z
        }

        let func_args = vec![1.0, 1.9, 1.0];
        let are_adjustings = vec![true, false, true];
        let derivative_funcs = derivative_funcs_numerical_difference_for_lm(&f, 2);
        let independent_vars = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let dependent_vars = vec![-0.9, 1.9, 7.3, 13.8, 23.5];
        let actual_vec = levenberg_marquardt(
            &f,
            func_args,
            are_adjustings,
            derivative_funcs,
            independent_vars,
            dependent_vars,
        );
        let tolerance = 1e-2;
        assert!((actual_vec[0] - (-0.945517241379294)).abs() < tolerance);
        assert_eq!(actual_vec[1], 1.9);
        assert!((actual_vec[2] - 1.04425287356321).abs() < tolerance);
    }
}
