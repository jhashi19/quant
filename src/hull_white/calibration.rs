use super::analysis;
use super::curve;
use super::optimization;

// pub fn swaption(
//     a: f64,
//     sigma: f64,
//     maturities: &Vec<f64>,
//     strike: f64,
//     swap_dates: &Vec<f64>,
//     op_type: analysis::SwaptionType,
//     curve: curve::Curve,
//     target_prices: &Vec<f64>,
// ) -> (f64, f64) {
//     let swap_dates = swap_dates.clone();
//     let op_type = op_type.clone();
//     let curve = curve.clone();
//     let strike = strike.clone();
//     let rap_swaption = move |maturity: f64, params: &[f64]| {
//         let [a, sigma] = [params[0], params[1]];
//         swaption_shifted_maturity(a, sigma, maturity, strike, &swap_dates, op_type, curve)
//     };
//     let func_args = vec![a, sigma];
//     let are_adjustings = vec![true, true];
//     let derivative_funcs =
//         optimization::derivative_funcs_numerical_difference_for_lm(&rap_swaption, 2);
//     let params = optimization::levenberg_marquardt(
//         &rap_swaption,
//         func_args,
//         are_adjustings,
//         derivative_funcs,
//         maturities.to_vec(),
//         target_prices.to_vec(),
//     );
//     (params[0], params[1])
// }

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
