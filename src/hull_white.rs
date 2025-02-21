mod analysis;
mod calibration;
mod curve;
mod data;
mod hw_lib;
mod interpolation;
mod math;
mod node;
mod optimization;
mod tree;

pub fn run() {
    /// Bermudan Swaptionの価格を計算する

    /// キャリブレーションのためのマーケットデータ
    // 満期
    let maturities = vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
    // 各満期共通のStrike
    let strikes = vec![0.005, 0.01, 0.015, 0.02, 0.025];
    // Swaptionの原資産のSwapのスケジュール
    let swap_dates = vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];
    // Swaptionの市場価格（満期×Strike）
    let market_prices = vec![
        vec![0.099, 0.090, 0.082, 0.078, 0.075, 0.070],
        vec![0.111, 0.100, 0.093, 0.089, 0.084, 0.081],
        vec![0.142, 0.139, 0.133, 0.128, 0.120, 0.118],
        vec![0.173, 0.169, 0.161, 0.158, 0.152, 0.149],
        vec![0.214, 0.210, 0.192, 0.190, 0.186, 0.181],
    ];

    // Hull-Whiteモデルのパラメータの初期値
    let init_a = 0.005;
    let init_sigma = 0.005;

    // キャリブレーション
    // maturitiesをfor文で回してHWのパラメータの値と区間の値をvecに追加していく
}
