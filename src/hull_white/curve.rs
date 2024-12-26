use super::data;
use super::interpolation::cubic_spline;
use super::optimization::{derivative_funcs_numerical_difference_for_lm, levenberg_marquardt};

#[derive(Clone, Copy, Debug)]
pub enum Curve {
    Ois,
    Libor6M,
    Libor12M,
}

pub fn match_curve(curve: Curve) -> (Vec<f64>, Vec<f64>) {
    match curve {
        Curve::Ois => (data::get_ois_date(), data::get_ois_rate()),
        Curve::Libor6M => (data::get_libor6m_date(), data::get_libor6m_rate()),
        Curve::Libor12M => (data::get_libor12m_date(), data::get_libor12m_rate()),
    }
}

pub fn df(curve: Curve, t: f64) -> f64 {
    let (dates, rates) = match_curve(curve);
    cubic_spline(&dates, &rates, t)
}

/// Svenssonによる瞬間フォワードレートの推定値を返します。
pub fn instantaneous_forward_rate(curve: Curve, t: f64) -> f64 {
    let (mut dates, mut rates) = match_curve(curve);
    // 先頭のデータは0日目のためDFからzero rateに変換するときに0割りが発生するため削除する。
    dates = dates[1..].to_vec();
    rates = rates[1..]
        .to_vec()
        .iter()
        .zip(dates.iter())
        .map(|(rate, date)| -rate.ln() / date)
        .collect();
    let beta0 = 0.1;
    let beta1 = 0.0;
    let beta2 = 0.01;
    let beta3 = 0.0;
    let tau1 = 2.0;
    let tau2 = 1.5;
    let params = vec![beta0, beta1, beta2, beta3, tau1, tau2];
    let are_adjustings = params.iter().map(|_| true).collect::<Vec<bool>>();
    let derivative_funcs =
        derivative_funcs_numerical_difference_for_lm(zero_rate_svensson, params.len());

    let adjusted_params = levenberg_marquardt(
        &zero_rate_svensson,
        params,
        are_adjustings,
        derivative_funcs,
        dates,
        rates,
    );

    instantaneous_forward_rate_svensson(t, &adjusted_params)
}

/// Svenssonによる瞬間フォワードレートの推定に使用するパラメータ算出のためのゼロレート関数
fn instantaneous_forward_rate_svensson(t: f64, params: &[f64]) -> f64 {
    let beta0 = params[0];
    let beta1 = params[1];
    let beta2 = params[2];
    let beta3 = params[3];
    let tau1 = params[4];
    let tau2 = params[5];
    beta0
        + beta1 * (-t / tau1).exp()
        + beta2 * (t / tau1) * (-t / tau1).exp()
        + beta3 * (t / tau2) * (-t / tau2).exp()
}

fn zero_rate_svensson(t: f64, params: &[f64]) -> f64 {
    let beta0 = params[0];
    let beta1 = params[1];
    let beta2 = params[2];
    let beta3 = params[3];
    let tau1 = params[4];
    let tau2 = params[5];
    beta0
        + beta1 * (1.0 - (-t / tau1).exp()) / (t / tau1)
        + beta2 * ((1.0 - (-t / tau1).exp()) / (t / tau1) - (-t / tau1).exp())
        + beta3 * ((1.0 - (-t / tau2).exp()) / (t / tau2) - (-t / tau2).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_instaneous_forward_rate() {
        let curve = Curve::Ois;
        let t = 1.0;
        let rate = instantaneous_forward_rate(curve, t);
        let tolerance = 1e-4;
        assert!((rate - 0.00581).abs() < tolerance);
    }
}
