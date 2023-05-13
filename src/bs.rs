use self::black_scholes::black_scholes;

mod black_scholes;

pub fn run() {
    let input = black_scholes::CalcInput {
        zero_rate: 0.02,
        vol: 0.2,
        term_annu: 0.5,
        strike: 98.0,
        underlying: 100.0,
    };

    let option_type = black_scholes::OptionType::Call;

    let opt_val = black_scholes(&input, option_type);
    println!("(analytical)price of european call option: {}", opt_val);
}
