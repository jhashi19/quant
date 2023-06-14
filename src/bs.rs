use self::black_scholes::black_scholes;

mod black_scholes;

pub fn run() {
    let input = black_scholes::CalcInput {
        underlying: 62.0,
        strike: 60.0,
        vol: 0.2,
        zero_rate: 0.1,
        term_annu: 5.0 / 12.0,
    };

    let option_type = black_scholes::OptionType::Call;

    let opt_val = black_scholes(&input, option_type);
    println!("(analytical)price of european call option: {}", opt_val);
}
