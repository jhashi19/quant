use self::black_scholes::black_scholes;

mod black_scholes;

pub fn run() {
    let input = black_scholes::CalcInput {
        zero_rate: 0.05,
        vol: 0.2,
        term_annu: 1.0,
        strike: 100.0,
        underlying: 100.0,
    };

    let option_type = black_scholes::OptionType::Put;

    let opt_val = black_scholes(&input, option_type);
    println!("(analytical)price of european put option: {}", opt_val);
}
