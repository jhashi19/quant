#[derive(Clone, Debug)]
pub struct Node {
    pub trans_to: Vec<Self>,   // 遷移先のノード
    pub trans_from: Vec<Self>, // 遷移元のノード
    pub rate: f64,             // 当ノードから次のノードまでの期間の金利
    pub rate_fluc_mean: f64,   // 当ノードの金利の変化幅の平均(期待値)
    pub rate_fluc_var: f64,    // 当ノードの金利の変化幅の分散
    pub p_up: f64,             // 遷移確率(上昇)
    pub p_mid: f64,            // 遷移確率(変化なし)
    pub p_down: f64,           // 遷移確率(下落)
    pub arrow_debreu: f64,     // Arrow Debreu price
}

impl Node {
    pub fn new() -> Self {
        Self {
            trans_to: vec![],
            trans_from: vec![],
            rate: 0.0,
            rate_fluc_mean: 0.0,
            rate_fluc_var: 0.0,
            p_up: 0.0,
            p_mid: 0.0,
            p_down: 0.0,
            arrow_debreu: 0.0,
        }
    }

    /// 金利の変化幅の期待値を返します。
    pub fn calc_fluc_mean(a: &f64, rate: &f64, time_interval: f64) -> f64 {
        ((-a + time_interval).exp() - 1.0) * rate
    }

    /// 金利の変化幅の分散を返します。
    pub fn calc_fluc_var(a: &f64, sigma: &f64, time_interval: f64) -> f64 {
        sigma.powi(2) * (1.0 - (-2.0 * a * time_interval).exp()) / (2.0 * a)
    }
}
