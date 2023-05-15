pub struct CalcInput {
    pub underlying: f64,
    pub strike: f64,
    pub vol: f64,
    pub zero_rate: f64,
    pub term_annu: f64,
}

pub fn crr_euro_call_layer(input: &CalcInput) {
    let CalcInput {
        underlying,
        strike,
        vol,
        zero_rate,
        term_annu,
    } = *input;

    let grid_num = 181;

    let delta_t = term_annu / grid_num as f64;
    let df = (-zero_rate * delta_t).exp();
    let (val_up, val_down) = crr(vol, delta_t);
    let rnp = ((zero_rate * delta_t).exp() - val_down) / (val_up - val_down);

    // create path
    let mut vals = vec![0.0; grid_num as usize];
    let val_up_square = val_up.powi(2);
    vals[0] = underlying * val_down.powi(grid_num - 1);
    for i in 1..grid_num as usize {
        vals[i] = vals[i - 1] * val_up_square;
    }

    // calc payoff
    let call_payoff = |x: &f64| (x - strike).max(0.0);
    let vals = vals.iter().map(call_payoff).collect::<Vec<f64>>();

    // backward induction
    let val = backward_euro(vals, grid_num, df, rnp);

    println!("(lattice layer)price of european call option: {}", val);

    // return (val_up, val_down)
    fn crr(vol: f64, delta_t: f64) -> (f64, f64) {
        let delta_t_sqrt = delta_t.sqrt();
        let val_up = (vol * delta_t_sqrt).exp();
        let val_down = (-vol * delta_t_sqrt).exp();
        (val_up, val_down)
    }

    fn backward_euro(mut vals: Vec<f64>, grid_num: i32, df: f64, rnp: f64) -> f64 {
        for i in (0..grid_num as usize).rev() {
            for j in 0..i {
                vals[j] = df * (rnp * vals[j + 1] + (1.0 - rnp) * vals[j]);
            }
        }
        vals[0]
    }
}
