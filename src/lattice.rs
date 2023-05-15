mod crr_adhoc;
mod crr_layer;

pub fn run() {
    crr_adhoc::crr_euro_call();

    let input = crr_layer::CalcInput {
        underlying: 100.0,
        strike: 98.0,
        vol: 0.2,
        zero_rate: 0.02,
        term_annu: 0.5,
    };
    crr_layer::crr_euro_call_layer(&input);
}
