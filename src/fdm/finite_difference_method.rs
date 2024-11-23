pub struct CalcInput {
    pub underlying: f64,
    pub strike: f64,
    pub vol: f64,
    pub zero_rate: f64,
    pub term_annu: f64,
}

// 有限差分法の陽解法でBlack-Scholes偏微分方程式の数値解を導出する。
pub fn explicit_fdm_bs(
    input: &CalcInput,
    p_max: f64,
    num_price_idx: usize,
    num_time_idx: usize,
) -> f64 {
    let CalcInput {
        underlying,
        strike,
        vol,
        zero_rate,
        term_annu,
    } = *input;

    let p_delta = p_max / num_price_idx as f64;
    let t_delta = term_annu / num_time_idx as f64;

    let mut opt_vals = vec![0.0; num_price_idx + 1];
    for idx in 0..num_price_idx + 1 {
        opt_vals[idx] = (p_delta * idx as f64 - strike).max(0.0);
    }

    // backwardで価格を計算し、もとのベクタを上書きしていく。
    for t_idx in (0..num_time_idx).rev() {
        // 複数要素から新たなベクタの要素を計算するため、一旦ベクタを複製してその要素を使って元のベクタを上書きする。
        let opt_vals_tmp = opt_vals.clone();
        for p_idx in 1..num_price_idx {
            let a = 0.5 * (zero_rate * p_idx as f64 + (vol * p_idx as f64).powi(2)) * t_delta;
            let b = 1.0 - ((vol * p_idx as f64).powi(2) + zero_rate) * t_delta;
            let c = 0.5 * (-zero_rate * p_idx as f64 + (vol * p_idx as f64).powi(2)) * t_delta;

            opt_vals[p_idx] =
                a * opt_vals_tmp[p_idx + 1] + b * opt_vals_tmp[p_idx] + c * opt_vals_tmp[p_idx - 1];
        }
        opt_vals[0] = 0.0;
        opt_vals[num_price_idx] =
            p_max - (-zero_rate * (term_annu - t_idx as f64 * t_delta)).exp() * strike;
    }

    // ひとまず簡単のためにunderlyingは整数で表せるものとする。
    // そのうち任意の実数に対応できるように修正する。(線形補間)
    opt_vals[(underlying * num_price_idx as f64 / p_max) as usize]
}

// Crank-Nicolson法でBlack-Scholes偏微分方程式の数値解を導出する。
pub fn crank_nicolson_fdm_bs(
    input: &CalcInput,
    p_max: f64,
    num_price_idx: usize,
    num_time_idx: usize,
) -> f64 {
    let CalcInput {
        underlying,
        strike,
        vol,
        zero_rate,
        term_annu,
    } = *input;

    let p_delta = p_max / num_price_idx as f64;
    let t_delta = term_annu / num_time_idx as f64;

    let vec_len = num_price_idx - 1;
    let mut p_vec: Vec<f64> = vec![0.0; vec_len];

    // backwardで価格を計算する際の初期値
    for j in 0..vec_len {
        p_vec[j] = ((j as f64 + 1.0) * p_delta - strike).max(0.0);
    }

    // Thomas法のためのセットアップ
    // メモリ節約のため、行列ではなく三行の対角成分のみのベクタで計算する。
    let mut upper_diag: Vec<f64> = vec![0.0; vec_len]; // 上側対角成分
    let mut middle_diag: Vec<f64> = vec![0.0; vec_len]; // 対角成分
    let mut lower_diag: Vec<f64> = vec![0.0; vec_len]; // 下側対角成分
    let mut rhs_coeff: Vec<f64> = vec![0.0; vec_len]; // CN法で連立方程式の右辺に現れる、対角成分に対応する係数

    for idx in 0..vec_len {
        upper_diag[idx] =
            0.25 * (-zero_rate * (idx + 1) as f64 - (vol * (idx + 1) as f64).powi(2)) * t_delta;
        middle_diag[idx] = 1.0 + 0.5 * ((vol * (idx + 1) as f64).powi(2) + zero_rate) * t_delta;
        lower_diag[idx] =
            0.25 * (zero_rate * (idx + 1) as f64 - (vol * (idx + 1) as f64).powi(2)) * t_delta;
        rhs_coeff[idx] = middle_diag[idx] - 2.0;
    }

    // 境界条件(価格の下限) 常に0とするためループの外側で定義している。
    let p_0 = 0.0;
    for i in (0..num_time_idx).rev() {
        // 境界条件(価格の上限)
        let p_vec_len = p_max - (-zero_rate * (term_annu - i as f64 * t_delta)).exp() * strike;
        let p_vec_len_next =
            p_max - (-zero_rate * (term_annu - (i as f64 + 1.0) * t_delta)).exp() * strike;

        // 前後のインデックスの値を使って計算するため、値が混ざらないようにp_vecの値をコピーして使用する。
        let p_vec_tmp = p_vec.clone();
        for j in 1..vec_len - 1 {
            p_vec[j] = -upper_diag[j] * p_vec_tmp[j + 1]
                - rhs_coeff[j] * p_vec_tmp[j]
                - lower_diag[j] * p_vec_tmp[j - 1];
        }
        p_vec[0] = -upper_diag[0] * p_vec_tmp[1]
            - rhs_coeff[0] * p_vec_tmp[0]
            - lower_diag[0] * p_0
            - lower_diag[0] * p_0;
        p_vec[vec_len - 1] = -upper_diag[vec_len - 1] * p_vec_len_next
            - rhs_coeff[vec_len - 1] * p_vec_tmp[vec_len - 1]
            - lower_diag[vec_len - 1] * p_vec_tmp[vec_len - 2]
            - upper_diag[vec_len - 1] * p_vec_len;

        p_vec = solve_by_thomas(&upper_diag, &middle_diag, &lower_diag, p_vec);
    }
    let mut price_vec = vec![0.0];
    price_vec.append(&mut p_vec);
    price_vec[(underlying * num_price_idx as f64 / p_max) as usize]
}

fn solve_by_thomas(
    upper: &Vec<f64>,
    middle: &Vec<f64>,
    lower: &Vec<f64>,
    mut p_vec: Vec<f64>,
) -> Vec<f64> {
    let mut upper = upper.clone();
    let mut middle = middle.clone();
    upper[0] /= middle[0];
    p_vec[0] /= middle[0];
    middle[0] = 1.0;
    for j in 1..upper.len() {
        upper[j] /= middle[j] - lower[j] * upper[j - 1];
        p_vec[j] = (p_vec[j] - p_vec[j - 1] * lower[j]) / (middle[j] - upper[j - 1] * lower[j]);
        middle[j] = 1.0;
    }

    for j in (0..upper.len() - 1).rev() {
        p_vec[j] -= upper[j] * p_vec[j + 1];
    }
    p_vec
}
