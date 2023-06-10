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
