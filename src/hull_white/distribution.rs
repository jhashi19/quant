use libm::erf;
use std::f64::{INFINITY, NEG_INFINITY};

pub fn std_normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / 2_f64.sqrt()))
}

// Acklam's algorithm
pub fn inverse_std_normal_cdf(x: f64) -> f64 {
    let a1 = -3.969683028665376e+01;
    let a2 = 2.209460984245205e+02;
    let a3 = -2.759285104469687e+02;
    let a4 = 1.383577518672690e+02;
    let a5 = -3.066479806614716e+01;
    let a6 = 2.506628277459239e+00;

    let b1 = -5.447609879822406e+01;
    let b2 = 1.615858368580409e+02;
    let b3 = -1.556989798598866e+02;
    let b4 = 6.680131188771972e+01;
    let b5 = -1.328068155288572e+01;

    let c1 = -7.784894002430293e-03;
    let c2 = -3.223964580411365e-01;
    let c3 = -2.400758277161838e+00;
    let c4 = -2.549732539343734e+00;
    let c5 = 4.374664141464968e+00;
    let c6 = 2.938163982698783e+00;

    let d1 = 7.784695709041462e-03;
    let d2 = 3.224671290700398e-01;
    let d3 = 2.445134137142996e+00;
    let d4 = 3.754408661907416e+00;

    // 計算式を分ける範囲の境界
    let x_low = 0.02425;
    let x_high = 1.0 - x_low;

    // 中央　中央の範囲に入ることが多いので、分岐の最初に持ってくることで判定処理を減らす
    if x_low <= x && x <= x_high {
        let q = x - 0.5;
        let r = q * q;
        return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q
            / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
    // 下側
    } else if 0.0 < x && x < x_low {
        let q = (-2.0 * x.ln()).sqrt();
        return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6)
            / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    // 上側
    } else if x_high < x && x < 1.0 {
        let q = (-2.0 * (1.0 - x).ln()).sqrt();
        return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6)
            / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    // 両端
    } else {
        if x >= 1.0 {
            return INFINITY;
        } else {
            return NEG_INFINITY;
        }
    }
}

pub fn wichura_inverse_normal_cdf(x: f64) -> f64 {
    let split1 = 0.425;
    let split2 = 5.0;
    let const1 = 0.180625;
    let const2 = 1.6;

    let a1 = 2.5090809287301226727e+03;
    let a2 = 3.3430575583588128105e+04;
    let a3 = 6.7265770927008700853e+04;
    let a4 = 4.5921953931549871457e+04;
    let a5 = 1.3731693765509461125e+04;
    let a6 = 1.9715909503065514427e+03;
    let a7 = 1.3314166789178437745e+02;
    let a8 = 3.3871328727963666080e+00;

    let b1 = 5.2264952788528545610e+03;
    let b2 = 2.8729085735721942674e+04;
    let b3 = 3.9307895800092710610e+04;
    let b4 = 2.1213794301586595867e+04;
    let b5 = 5.3941960214247511077e+03;
    let b6 = 6.8718700749205790830e+02;
    let b7 = 4.2313330701600911252e+01;

    let c1 = 7.74545014278341407640e-04;
    let c2 = 2.27238449892691845833e-02;
    let c3 = 2.41780725177450611770e-01;
    let c4 = 1.27045825245236838258e+00;
    let c5 = 3.64784832476320460504e+00;
    let c6 = 5.76949722146069140550e+00;
    let c7 = 4.63033784615654529590e+00;
    let c8 = 1.42343711074968357734e+00;

    let d1 = 1.05075007164441684324e-09;
    let d2 = 5.47593808499534494600e-04;
    let d3 = 1.51986665636164571966e-02;
    let d4 = 1.48103976427480074590e-01;
    let d5 = 6.89767334985100004550e-01;
    let d6 = 1.67638483018380384940e+00;
    let d7 = 2.05319162663775882187e+00;

    let e1 = 2.01033439929228813265e-07;
    let e2 = 2.71155556874348757815e-05;
    let e3 = 1.24266094738807843860e-03;
    let e4 = 2.65321895265761230930e-02;
    let e5 = 2.96560571828504891230e-01;
    let e6 = 1.78482653991729133580e+00;
    let e7 = 5.46378491116411436990e+00;
    let e8 = 6.65790464350110377720e+00;

    let f1 = 2.04426310338993978564e-15;
    let f2 = 1.42151175831644588870e-07;
    let f3 = 1.84631831751005468180e-05;
    let f4 = 7.86869131145613259100e-04;
    let f5 = 1.48753612908506148525e-02;
    let f6 = 1.36929880922735805310e-01;
    let f7 = 5.99832206555887937690e-01;

    let mut r;

    let q = x - 0.5;
    if q.abs() <= split1 {
        let r = const1 - q * q;
        return q * (((((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * r + a7) * r + a8)
            / (((((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + b6) * r + b7) * r + 1.0);
    } else {
        if q < 0.0 {
            r = x;
        } else {
            r = 1.0 - x;
        }
        if r <= 0.0 {
            return 0.0;
        }
        r = (-r.ln()).sqrt();
        if r <= split2 {
            r -= const2;
            return q.signum()
                * (((((((c1 * r + c2) * r + c3) * r + c4) * r + c5) * r + c6) * r + c7) * r + c8)
                / (((((((d1 * r + d2) * r + d3) * r + d4) * r + d5) * r + d6) * r + d7) * r + 1.0);
        } else {
            r -= split2;
            return q.signum()
                * (((((((e1 * r + e2) * r + e3) * r + e4) * r + e5) * r + e6) * r + e7) * r + e8)
                / (((((((f1 * r + f2) * r + f3) * r + f4) * r + f5) * r + f6) * r + f7) * r + 1.0);
        }
    }
}

pub fn moro_inverse_normal_cdf(average: f64, std_dev: f64, x: f64) -> f64 {
    let a0: f64 = 2.50662823884;
    let a1: f64 = -18.61500062529;
    let a2: f64 = 41.39119773534;
    let a3: f64 = -25.44106049637;

    let b0: f64 = -8.47351093090;
    let b1: f64 = 23.08336743743;
    let b2: f64 = -21.06224101826;
    let b3: f64 = 3.13082909833;

    let c0: f64 = 0.3374754822726147;
    let c1: f64 = 0.9761690190917186;
    let c2: f64 = 0.1607979714918209;
    let c3: f64 = 0.0276438810333863;
    let c4: f64 = 0.0038405729373609;
    let c5: f64 = 0.0003951896511919;
    let c6: f64 = 0.0000321767881768;
    let c7: f64 = 0.0000002888167364;
    let c8: f64 = 0.0000003960315187;

    let mut result: f64;
    let tmp = x - 0.5;

    if tmp.abs() < 0.42 {
        result = tmp.powi(2);
        result = tmp * (((a3 * result + a2) * result + a1) * result + a0)
            / ((((b3 * result + b2) * result + b1) * result + b0) * result + 1.0);
    } else {
        if x < 0.5 {
            result = x;
        } else {
            result = 1.0 - x;
        }
        result = (-result.ln()).ln();
        result = c0
            + result
                * (c1
                    + result
                        * (c2
                            + result
                                * (c3
                                    + result
                                        * (c4
                                            + result
                                                * (c5
                                                    + result
                                                        * (c6 + result * (c7 + result * c8)))))));
        if x < 0.5 {
            result = -result;
        }
    }
    average + result * std_dev
}

pub fn moro_inverse_std_normal_cdf(x: f64) -> f64 {
    moro_inverse_normal_cdf(0.0, 1.0, x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wichura_inverse_normal_cdf() {
        let threshold = 10_f64.powi(-13);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.5398278372770290) - (-0.1)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.5792597094391030) - (-0.2)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.6179114221889526) - (-0.3)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.6554217416103242) - (-0.4)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.6914624612740131) - (-0.5)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.7257468822499270) - (-0.6)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.7580363477769270) - (-0.7)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.7881446014166033) - (-0.8)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.8159398746532405) - (-0.9)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.8413447460685429) - (-1.0)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.9331927987311419) - (-1.5)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.9772498680518208) - (-2.0)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.9937903346742240) - (-2.5)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.9986501019683699) - (-3.0)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.9997673709209645) - (-3.5)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(1.0 - 0.9999683287581669) - (-4.0)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.5)).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.5398278372770290) - 0.1).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.5792597094391030) - 0.2).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.6179114221889526) - 0.3).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.6554217416103242) - 0.4).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.6914624612740131) - 0.5).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.7257468822499270) - 0.6).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.7580363477769270) - 0.7).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.7881446014166033) - 0.8).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.8159398746532405) - 0.9).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.8413447460685429) - 1.0).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.9331927987311419) - 1.5).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.9772498680518208) - 2.0).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.9937903346742240) - 2.5).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.9986501019683699) - 3.0).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.9997673709209645) - 3.5).abs() < threshold);
        assert!((wichura_inverse_normal_cdf(0.9999683287581669) - 4.0).abs() < threshold);
    }

    #[test]
    fn test_inverse_std_normal_cdf() {
        let threshold = 10_f64.powi(-8);
        assert!((inverse_std_normal_cdf(1.0 - 0.5398278372770290) - (-0.1)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.5792597094391030) - (-0.2)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.6179114221889526) - (-0.3)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.6554217416103242) - (-0.4)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.6914624612740131) - (-0.5)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.7257468822499270) - (-0.6)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.7580363477769270) - (-0.7)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.7881446014166033) - (-0.8)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.8159398746532405) - (-0.9)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.8413447460685429) - (-1.0)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.9331927987311419) - (-1.5)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.9772498680518208) - (-2.0)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.9937903346742240) - (-2.5)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.9986501019683699) - (-3.0)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.9997673709209645) - (-3.5)).abs() < threshold);
        assert!((inverse_std_normal_cdf(1.0 - 0.9999683287581669) - (-4.0)).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.5)).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.5398278372770290) - 0.1).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.5792597094391030) - 0.2).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.6179114221889526) - 0.3).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.6554217416103242) - 0.4).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.6914624612740131) - 0.5).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.7257468822499270) - 0.6).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.7580363477769270) - 0.7).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.7881446014166033) - 0.8).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.8159398746532405) - 0.9).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.8413447460685429) - 1.0).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.9331927987311419) - 1.5).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.9772498680518208) - 2.0).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.9937903346742240) - 2.5).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.9986501019683699) - 3.0).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.9997673709209645) - 3.5).abs() < threshold);
        assert!((inverse_std_normal_cdf(0.9999683287581669) - 4.0).abs() < threshold);
    }

    #[test]
    fn test_moro_inverse_std_normal_cdf() {
        let threshold = 10_f64.powi(-8);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.5398278372770290) - (-0.1)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.5792597094391030) - (-0.2)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.6179114221889526) - (-0.3)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.6554217416103242) - (-0.4)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.6914624612740131) - (-0.5)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.7257468822499270) - (-0.6)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.7580363477769270) - (-0.7)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.7881446014166033) - (-0.8)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.8159398746532405) - (-0.9)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.8413447460685429) - (-1.0)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.9331927987311419) - (-1.5)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.9772498680518208) - (-2.0)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.9937903346742240) - (-2.5)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.9986501019683699) - (-3.0)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.9997673709209645) - (-3.5)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(1.0 - 0.9999683287581669) - (-4.0)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.5)).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.5398278372770290) - 0.1).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.5792597094391030) - 0.2).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.6179114221889526) - 0.3).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.6554217416103242) - 0.4).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.6914624612740131) - 0.5).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.7257468822499270) - 0.6).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.7580363477769270) - 0.7).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.7881446014166033) - 0.8).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.8159398746532405) - 0.9).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.8413447460685429) - 1.0).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.9331927987311419) - 1.5).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.9772498680518208) - 2.0).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.9937903346742240) - 2.5).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.9986501019683699) - 3.0).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.9997673709209645) - 3.5).abs() < threshold);
        assert!((moro_inverse_std_normal_cdf(0.9999683287581669) - 4.0).abs() < threshold);
    }

    #[test]
    fn test_std_normal_cdf() {
        let threshold = 10_f64.powi(-15);
        assert!((std_normal_cdf(-0.1) - (1.0 - 0.5398278372770290)).abs() < threshold);
        assert!((std_normal_cdf(-0.2) - (1.0 - 0.5792597094391030)).abs() < threshold);
        assert!((std_normal_cdf(-0.3) - (1.0 - 0.6179114221889526)).abs() < threshold);
        assert!((std_normal_cdf(-0.4) - (1.0 - 0.6554217416103242)).abs() < threshold);
        assert!((std_normal_cdf(-0.5) - (1.0 - 0.6914624612740131)).abs() < threshold);
        assert!((std_normal_cdf(-0.6) - (1.0 - 0.7257468822499270)).abs() < threshold);
        assert!((std_normal_cdf(-0.7) - (1.0 - 0.7580363477769270)).abs() < threshold);
        assert!((std_normal_cdf(-0.8) - (1.0 - 0.7881446014166033)).abs() < threshold);
        assert!((std_normal_cdf(-0.9) - (1.0 - 0.8159398746532405)).abs() < threshold);
        assert!((std_normal_cdf(-1.0) - (1.0 - 0.8413447460685429)).abs() < threshold);
        assert!((std_normal_cdf(-1.5) - (1.0 - 0.9331927987311419)).abs() < threshold);
        assert!((std_normal_cdf(-2.0) - (1.0 - 0.9772498680518208)).abs() < threshold);
        assert!((std_normal_cdf(-2.5) - (1.0 - 0.9937903346742240)).abs() < threshold);
        assert!((std_normal_cdf(-3.0) - (1.0 - 0.9986501019683699)).abs() < threshold);
        assert!((std_normal_cdf(-3.5) - (1.0 - 0.9997673709209645)).abs() < threshold);
        assert!((std_normal_cdf(-4.0) - (1.0 - 0.9999683287581669)).abs() < threshold);
        assert!((std_normal_cdf(0.0) - 0.5).abs() < threshold);
        assert!((std_normal_cdf(0.1) - 0.5398278372770290).abs() < threshold);
        assert!((std_normal_cdf(0.2) - 0.5792597094391030).abs() < threshold);
        assert!((std_normal_cdf(0.3) - 0.6179114221889526).abs() < threshold);
        assert!((std_normal_cdf(0.4) - 0.6554217416103242).abs() < threshold);
        assert!((std_normal_cdf(0.5) - 0.6914624612740131).abs() < threshold);
        assert!((std_normal_cdf(0.6) - 0.7257468822499270).abs() < threshold);
        assert!((std_normal_cdf(0.7) - 0.7580363477769270).abs() < threshold);
        assert!((std_normal_cdf(0.8) - 0.7881446014166033).abs() < threshold);
        assert!((std_normal_cdf(0.9) - 0.8159398746532405).abs() < threshold);
        assert!((std_normal_cdf(1.0) - 0.8413447460685429).abs() < threshold);
        assert!((std_normal_cdf(1.5) - 0.9331927987311419).abs() < threshold);
        assert!((std_normal_cdf(2.0) - 0.9772498680518208).abs() < threshold);
        assert!((std_normal_cdf(2.5) - 0.9937903346742240).abs() < threshold);
        assert!((std_normal_cdf(3.0) - 0.9986501019683699).abs() < threshold);
        assert!((std_normal_cdf(3.5) - 0.9997673709209645).abs() < threshold);
        assert!((std_normal_cdf(4.0) - 0.9999683287581669).abs() < threshold);
    }
}
