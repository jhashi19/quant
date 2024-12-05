use super::curve::{
    self,
    Curve::{self, Ois},
};
use super::math::std_normal_cdf;

/* One Factor Hull White
dr_t = (θ_t - a * r_t) * dt - σ * dW_t */

// Day Count Conventionは考慮しない。
#[derive(Clone, Copy, Debug)]
pub enum OptionType {
    Call,
    Put,
}

#[derive(Clone, Copy, Debug)]
pub enum CapFloorType {
    Cap,
    Floor,
}

/// 0時点までのディスカウントファクターを返します。
/// * `curve` - ディスカウントカーブの種類(OIS、LIBOR6M、LIBOR12M)
/// * `t` - 割り戻すCFの時点
pub fn df(curve: curve::Curve, t: f64) -> f64 {
    curve::df(curve, t)
}

/// 割引債オプションの理論価格を返します。
/// * `mat_u` - 原資産の債券の満期日
/// * `mat_o` - オプションの満期日
/// * `strike` - 権利行使価格
/// * `vol` - ボラティリティ
/// * `op_type` - Call/Put
pub fn discount_bond_option(
    mat_u: f64,
    mat_o: f64,
    strike: f64,
    vol: f64,
    op_type: OptionType,
) -> f64 {
    let sign = match op_type {
        OptionType::Call => 1.0,
        OptionType::Put => -1.0,
    };
    let d = (df(Ois, mat_u) / (strike * df(Ois, mat_o))).ln() + 0.5 * vol.powi(2) / vol;
    sign * df(Ois, mat_u) * std_normal_cdf(sign * d)
        - sign * df(Ois, mat_o) * std_normal_cdf(sign * (d - vol))
}

/// 割引債オプションのボラティリティを返します。
pub fn dbo_sigma(a: &f64, sigma: &f64, mat_u: f64, mat_o: f64) -> f64 {
    sigma / a
        * (1.0 - (-2.0 * a * mat_o).exp() / (2.0 * a)).powf(0.5)
        * (1.0 - (-a * (mat_u - mat_o)).exp())
}

/// Caplet、Floorletの理論価格を返します。
pub fn capfloorlet(
    date_s: f64,
    date_e: f64,
    strike: f64,
    vol: f64,
    cf_type: CapFloorType,
    curve: Curve,
) -> f64 {
    let sign = match cf_type {
        CapFloorType::Cap => -1.0,
        CapFloorType::Floor => 1.0,
    };
    let delta = date_e - date_s;
    let d = (((1.0 + delta * strike) * df(curve, date_e) / df(curve, date_s)).ln()
        + 0.5 * vol.powi(2))
        / vol;
    df(Ois, date_e)
        * (sign * (1.0 + delta * strike) * std_normal_cdf(sign * d)
            - sign * df(curve, date_s) / df(curve, date_e) * std_normal_cdf(sign * (d - vol)))
}

/// Caplet、Floorletのボラティリティを返します。
pub fn capfloorlet_sigma(a: &f64, sigma: &f64, date_s: f64, date_e: f64) -> f64 {
    sigma / a * (1.0 - (-2.0 * date_s).exp()) / (0.5 * a).powf(0.5)
        * (1.0 - (-a * (date_e - date_s)))
}

/// Cap、Floorの理論価格を返します。
pub fn capfloor(
    dates: &Vec<f64>, // 各Caplet/Floorletの参照レートのスタートとエンドの日付(エンドが次のCFのスタートと一致すると仮定)
    vols: &Vec<f64>,  // 各Caplet/Floorletのボラティリティ
    strike: f64,
    cf_type: CapFloorType,
    curve: Curve,
) -> f64 {
    if dates.len() != vols.len() {
        panic!("CFの数とボラティリティの数が合っていません");
    }
    let mut price = 0.0;
    for i in 0..dates.len() {
        price += capfloorlet(dates[i], dates[i + 1], strike, vols[i], cf_type, curve);
    }
    price
}

/// Swaptionの理論価格を返します。
pub fn swaption() -> f64 {
    0.0
}
