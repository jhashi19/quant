pub fn crr_euro_call() {
    // set parameter
    let init_und = 100.0;
    let strike = 98.0;
    let vol = 0.2;
    let zero_rate = 0.02;
    let mat = 0.5;
    let grid_num = 181;

    // calc parameter
    let delta_t = mat / grid_num as f64;
    let delta_t_sqrt = delta_t.sqrt();
    let val_up = (vol * delta_t_sqrt).exp();
    let val_down = (-vol * delta_t_sqrt).exp();
    let df = (-zero_rate * delta_t).exp();
    let rnp = ((zero_rate * delta_t).exp() - val_down) / (val_up - val_down);

    // create path
    let mut vals = vec![0.0; grid_num as usize];

    let val_up_square = val_up.powi(2);
    vals[0] = init_und * val_down.powi(grid_num - 1);
    for i in 1..grid_num as usize {
        vals[i] = vals[i - 1] * val_up_square;
    }

    let mut vals = vals
        .iter()
        .map(|x| (x - strike).max(0.0))
        .collect::<Vec<f64>>();

    for i in (0..grid_num as usize).rev() {
        for j in 0..i {
            vals[j] = df * (rnp * vals[j + 1] + (1.0 - rnp) * vals[j]);
        }
    }
    println!("underlying price: {}", init_und);
    println!("strike price: {}", strike);
    println!("volatility: {}", vol);
    println!("zero rate: {}", zero_rate);
    println!("maturity: {}", mat);
    println!("gird number: {}", grid_num);
    println!("price of european call option: {}", vals[0]);
}
