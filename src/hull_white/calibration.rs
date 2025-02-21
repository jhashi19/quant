use super::analysis;
use super::curve;
use super::optimization;

pub fn swaption_with_maturities(
    init_a: f64,
    init_sigma: f64,
    maturities: &Vec<f64>,
    strike: f64,
    swap_dates: &Vec<f64>,
    op_type: analysis::SwaptionType,
    curve: curve::Curve,
    target_prices: &Vec<f64>,
) -> (f64, f64) {
    let rap_swaption = move |maturity: f64, params: &[f64]| {
        let [a, sigma] = [params[0], params[1]];
        swaption_shifted_maturity(a, sigma, maturity, strike, &swap_dates, op_type, curve)
    };
    let func_args = vec![init_a, init_sigma];
    let are_adjustings = vec![true, true];
    let derivative_funcs =
        optimization::derivative_funcs_numerical_difference_for_lm(&rap_swaption, 2);
    let adjusted_params = optimization::levenberg_marquardt(
        &rap_swaption,
        func_args,
        are_adjustings,
        derivative_funcs,
        maturities.to_vec(),
        target_prices.to_vec(),
    );
    (adjusted_params[0], adjusted_params[1])
}

fn swaption_shifted_maturity(
    a: f64,
    sigma: f64,
    mat_op: f64,
    strike: f64,
    swap_dates: &Vec<f64>,
    op_type: analysis::SwaptionType,
    curve: curve::Curve,
) -> f64 {
    let sliced_swap_dates = slice_swap_dates(mat_op, swap_dates);
    analysis::swaption(a, sigma, mat_op, strike, &sliced_swap_dates, op_type, curve)
}

fn slice_swap_dates(maturity: f64, swap_dates: &Vec<f64>) -> Vec<f64> {
    swap_dates
        .iter()
        .filter(|&date| date > &maturity)
        .map(|&x| x)
        .collect()
}

pub fn swaption_with_strikes(
    init_a: f64,
    init_sigma: f64,
    maturity: f64,
    strikes: &Vec<f64>,
    swap_dates: &Vec<f64>,
    op_type: analysis::SwaptionType,
    curve: curve::Curve,
    target_prices: &Vec<f64>,
) -> (f64, f64) {
    let rap_swaption = move |strike: f64, params: &[f64]| {
        let [a, sigma] = [params[0], params[1]];
        analysis::swaption(a, sigma, maturity, strike, &swap_dates, op_type, curve)
    };
    let func_args = vec![init_a, init_sigma];
    let are_adjustings = vec![true, true];
    let derivative_funcs =
        optimization::derivative_funcs_numerical_difference_for_lm(&rap_swaption, 2);
    let adjusted_params = optimization::levenberg_marquardt(
        &rap_swaption,
        func_args,
        are_adjustings,
        derivative_funcs,
        strikes.to_vec(),
        target_prices.to_vec(),
    );
    (adjusted_params[0], adjusted_params[1])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_slice_swap_dates() {
        let swap_dates = vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
        assert_eq!(
            slice_swap_dates(0.3, &swap_dates),
            vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        );
        assert_eq!(
            slice_swap_dates(0.49, &swap_dates),
            vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        );
        assert_eq!(
            slice_swap_dates(0.5, &swap_dates),
            vec![1.0, 1.5, 2.0, 2.5, 3.0]
        );
        assert_eq!(
            slice_swap_dates(0.51, &swap_dates),
            vec![1.0, 1.5, 2.0, 2.5, 3.0]
        );
        assert_eq!(slice_swap_dates(1.0, &swap_dates), vec![1.5, 2.0, 2.5, 3.0]);
        assert_eq!(slice_swap_dates(1.5, &swap_dates), vec![2.0, 2.5, 3.0]);
    }

    #[test]
    fn test_swaption_with_maturities() {
        let a = 0.005;
        let sigma = 0.005;
        let maturities = vec![0.5, 1.0, 1.5];
        let strike = 0.004;
        let swap_dates = vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];
        let op_type = analysis::SwaptionType::Payer;
        let curve = curve::Curve::Ois;
        let target_prices = vec![0.12, 0.14, 0.15];
        let (adjusted_a, adjusted_sigma) = swaption_with_maturities(
            a,
            sigma,
            &maturities,
            strike,
            &swap_dates,
            op_type,
            curve,
            &target_prices,
        );
        println!(
            "maturities[0]: {}",
            swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturities[0],
                strike,
                &swap_dates,
                op_type,
                curve
            )
        );
        assert!(
            (swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturities[0],
                strike,
                &swap_dates,
                op_type,
                curve
            ) - target_prices[0])
                .abs()
                < 1e-2
        );
        println!(
            "maturities[1]: {}",
            swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturities[1],
                strike,
                &swap_dates,
                op_type,
                curve
            )
        );
        assert!(
            (swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturities[1],
                strike,
                &swap_dates,
                op_type,
                curve
            ) - target_prices[1])
                .abs()
                < 1e-2
        );
        println!(
            "maturities[2]: {}",
            swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturities[2],
                strike,
                &swap_dates,
                op_type,
                curve
            )
        );
        assert!(
            (swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturities[2],
                strike,
                &swap_dates,
                op_type,
                curve
            ) - target_prices[2])
                .abs()
                < 1e-2
        );
    }

    #[test]
    fn test_swaption_with_strikes() {
        let a = 0.005;
        let sigma = 0.005;
        let maturity = 1.0;
        let strikes = vec![0.004, 0.007, 0.01];
        let swap_dates = vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];
        let op_type = analysis::SwaptionType::Payer;
        let curve = curve::Curve::Ois;
        let target_prices = vec![0.14, 0.13, 0.13];
        let (adjusted_a, adjusted_sigma) = swaption_with_strikes(
            a,
            sigma,
            maturity,
            &strikes,
            &swap_dates,
            op_type,
            curve,
            &target_prices,
        );
        println!(
            "maturities[0]: {}",
            swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturity,
                strikes[0],
                &swap_dates,
                op_type,
                curve
            )
        );
        assert!(
            (swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturity,
                strikes[0],
                &swap_dates,
                op_type,
                curve
            ) - target_prices[0])
                .abs()
                < 1e-2
        );
        println!(
            "maturities[1]: {}",
            swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturity,
                strikes[1],
                &swap_dates,
                op_type,
                curve
            )
        );
        assert!(
            (swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturity,
                strikes[1],
                &swap_dates,
                op_type,
                curve
            ) - target_prices[1])
                .abs()
                < 1e-2
        );
        println!(
            "maturities[2]: {}",
            swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturity,
                strikes[2],
                &swap_dates,
                op_type,
                curve
            )
        );
        assert!(
            (swaption_shifted_maturity(
                adjusted_a,
                adjusted_sigma,
                maturity,
                strikes[2],
                &swap_dates,
                op_type,
                curve
            ) - target_prices[2])
                .abs()
                < 1e-2
        );
    }
}
