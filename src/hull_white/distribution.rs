pub fn moro_inverse_cdf_normal(average: f64, std_dev: f64, x: f64) -> f64 {
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

pub fn moro_inverse_cdf_std_normal(x: f64) -> f64 {
    moro_inverse_cdf_normal(0.0, 1.0, x)
}
