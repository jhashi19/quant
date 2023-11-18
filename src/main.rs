use std::env;

mod bs;
mod fdm;
mod hull_white;
mod lattice;
mod lsm;
mod mc;
mod sabr;

fn main() {
    let args: Vec<String> = env::args().collect();
    let module = args[1].as_str();
    match module {
        "bs" => bs::run(),
        "fdm" => fdm::run(),
        "hull_white" => hull_white::run(),
        "lattice" => lattice::run(),
        "lsm" => lsm::run(),
        "mc" => mc::run(),
        "sabr" => sabr::run(),
        _ => println!("there is no module: {}", module),
    }
}
