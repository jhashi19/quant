use super::data;
use super::interpolation::cubic_spline;

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
