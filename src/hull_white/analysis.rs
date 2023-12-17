/* One Factor Hull White
dr_t = (θ_t - a * r_t) * dt - σ * dW_t */

use super::distribution::std_normal_cdf;

pub enum OptionType {
    Call,
    Put,
}

/// ディスカウントファクターを返します。
/// * `t` -  割戻先
/// * `matu` - 割戻元
pub fn df(t: &f64, matu: &f64) -> f64 {
    0.0
}

/// 割引債オプションの理論価格を返します。
/// * `a` - HWモデルのa
/// * `sigma` - HWモデルのsigma
/// * `t` - 計算基準日
/// * `mat_u` - 原資産の債券の満期日
/// * `mat_o` - オプションの満期日
/// * `strike` - 権利行使価格
/// * `op_type` - Call/Put
pub fn discount_bond_option(
    a: &f64,
    sigma: &f64,
    t: &f64,
    mat_u: &f64,
    mat_o: &f64,
    strike: &f64,
    op_type: OptionType,
) -> f64 {
    let cp_sign = match op_type {
        OptionType::Call => 1.0,
        OptionType::Put => -1.0,
    };
    let sigma_p = sigma / a
        * (1.0 - (-2.0 * a * (mat_u - t)).exp() / (2.0 * a)).powf(0.5)
        * (1.0 - (-a * (mat_o - mat_u)).exp());
    let d = (df(t, mat_u) / (strike * df(t, mat_o))).ln() + 0.5 * sigma_p.powi(2) / sigma_p;
    cp_sign * df(t, mat_u) * std_normal_cdf(cp_sign * d)
        - cp_sign * df(t, mat_o) * std_normal_cdf(cp_sign * (d - sigma_p))
}

/// Caplet、Floorletの理論価格を返します。
pub fn capfloorlet() -> f64 {
    // discount_bond_option使って計算する？(Hull-Whiteのpaper)
    0.0
}

/// Cap、Floorの理論価格を返します。
pub fn capfloor() -> f64 {
    0.0
}

/// Swaptionの理論価格を返します。
pub fn swaption() -> f64 {
    0.0
}
