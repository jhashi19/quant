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

#[derive(Clone, Copy, Debug)]
pub enum SwaptionType {
    Payer,
    Receiver,
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
pub fn dbo_vol(a: f64, sigma: f64, mat_u: f64, mat_o: f64) -> f64 {
    sigma / a
        * (1.0 - (-2.0 * a * mat_o).exp() / (2.0 * a)).sqrt()
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
pub fn capfloorlet_vol(a: f64, sigma: f64, date_s: f64, date_e: f64) -> f64 {
    sigma / a * (1.0 - (-2.0 * date_s).exp()) / (0.5 * a).sqrt() * (1.0 - (-a * (date_e - date_s)))
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
/// ATMとなるショートレートの逆算を省略するための近似値を返します。
pub fn swaption(
    swap_dates: &Vec<f64>,
    strike: f64,
    vol: f64,
    op_type: SwaptionType,
    curve: Curve,
) -> f64 {
    let sign = match op_type {
        SwaptionType::Payer => -1.0,
        SwaptionType::Receiver => 1.0,
    };

    let coupon_bearing_bond = coupon_bearing_bond(&swap_dates, strike, curve);

    // Swaption価格の計算
    let d = ((coupon_bearing_bond * df(curve, swap_dates[1])
        / (df(Ois, swap_dates[1]) * df(curve, swap_dates[0])))
    .ln()
        + 0.5 * vol.powi(2))
        / vol;
    sign * coupon_bearing_bond * std_normal_cdf(sign * d)
        - sign * df(Ois, swap_dates[1]) * df(curve, swap_dates[0]) / df(curve, swap_dates[1])
            * std_normal_cdf(sign * (d - vol))
}

pub fn swaption_vol(
    a: f64,
    sigma: f64,
    mat_op: f64,
    swap_dates: &Vec<f64>,
    strike: f64,
    curve: Curve,
) -> f64 {
    let b = |mat: f64| 1.0 / a * (1.0 - (-a * mat).exp());
    let coupons = coupons(swap_dates, strike, curve);
    let coupon_bearing_bond = coupon_bearing_bond(swap_dates, strike, curve);
    let weighted_b = coupons.iter().enumerate().fold(0.0, |acc, (i, coupon)| {
        acc + coupon * df(Ois, swap_dates[i + 1]) * (b(swap_dates[0]) - b(swap_dates[i + 1]))
    });
    let modified_mat =
        -1.0 / a * (1.0 - a * (b(swap_dates[0]) - 1.0 / coupon_bearing_bond * weighted_b)).ln();
    sigma / a
        * ((1.0 - (-2.0 * a * mat_op).exp()) / (2.0 * a)).sqrt()
        * ((-a * (swap_dates[0] - mat_op)).exp() - (-a * (modified_mat - mat_op)).exp())
}

fn coupons(swap_dates: &Vec<f64>, strike: f64, curve: Curve) -> Vec<f64> {
    let swap_length = swap_dates.len();
    let mut coupons: Vec<f64> = Vec::with_capacity(swap_length - 1);
    for i in 0..swap_length - 2 {
        coupons.push(
            1.0 - (df(curve, swap_dates[i + 1]) * df(Ois, swap_dates[i + 2]))
                / (df(curve, swap_dates[i + 2]) * df(Ois, swap_dates[i + 1]))
                + strike * (swap_dates[i + 1] - swap_dates[i]),
        );
    }
    coupons.push(1.0 + strike * (swap_dates[swap_length - 1] - swap_dates[swap_length - 2]));
    coupons
}

fn coupon_bearing_bond(swap_dates: &Vec<f64>, strike: f64, curve: Curve) -> f64 {
    coupons(swap_dates, strike, curve)
        .iter()
        .zip(swap_dates[1..].iter())
        .fold(0.0, |acc, (coupon, date)| acc + coupon * df(Ois, *date))
}
