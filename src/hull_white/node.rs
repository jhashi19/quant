#[derive(Clone, Debug)]
pub struct Node {
    pub transition_to: Vec<usize>, // 遷移先のノードの金利方向のセンタリングされたインデックスのベクタ(要素数は3)
    pub rate: f64,                 // 当ノードから次のノードまでの期間の金利
    pub rate_fluctuation_mean: f64, // 当ノードの金利の変化幅の平均(期待値)
    pub rate_fluctuation_var: f64, // 当ノードの金利の変化幅の分散
    pub prob_up: f64,              // 遷移確率(上昇)
    pub prob_mid: f64,             // 遷移確率(中間)
    pub prob_down: f64,            // 遷移確率(下落)
    pub arrow_debreu: f64,         // Arrow Debreu price
}

impl Node {
    pub fn new() -> Self {
        Self {
            transition_to: vec![],
            rate: 0.0,
            rate_fluctuation_mean: 0.0,
            rate_fluctuation_var: 0.0,
            prob_up: 0.0,
            prob_mid: 0.0,
            prob_down: 0.0,
            arrow_debreu: 0.0,
        }
    }

    /// 金利の変化幅の期待値を返します。
    pub fn calc_fluctuation_mean(a: &f64, rate: f64, time_interval: f64) -> f64 {
        ((-a + time_interval).exp() - 1.0) * rate
    }

    /// 金利の変化幅の分散を返します。
    pub fn calc_fluctuation_var(a: &f64, sigma: &f64, time_interval: f64) -> f64 {
        sigma.powi(2) * (1.0 - (-2.0 * a * time_interval).exp()) / (2.0 * a)
    }

    /// 遷移確率をtupleで返します。<br>
    /// 戻り値：(up, mid, down)
    pub fn calc_transition_prob(
        rate: f64,
        rate_fluctuation_mean: f64,
        rate_fluctuation_var: f64,
        next_node_mid_rate: f64,
        next_rate_fluctuation: f64,
    ) -> (f64, f64, f64) {
        let alpha = Self::calc_alpha(
            rate,
            rate_fluctuation_mean,
            next_node_mid_rate,
            next_rate_fluctuation,
        );
        (
            Self::calc_prob_up(alpha, rate_fluctuation_var, next_rate_fluctuation),
            Self::calc_prob_mid(alpha, rate_fluctuation_var, next_rate_fluctuation),
            Self::calc_prob_down(alpha, rate_fluctuation_var, next_rate_fluctuation),
        )
    }

    /// 各遷移確率の計算に使用する値αを返します。
    fn calc_alpha(
        rate: f64,
        rate_fluctuation_mean: f64,
        next_node_mid_rate: f64,
        next_rate_fluctuation: f64,
    ) -> f64 {
        (rate + rate_fluctuation_mean - next_node_mid_rate) / next_rate_fluctuation
    }

    /// 遷移確率(上昇)を返します。
    fn calc_prob_up(alpha: f64, rate_fluctuation_var: f64, next_rate_fluctuation: f64) -> f64 {
        rate_fluctuation_var / (2.0 * next_rate_fluctuation.powi(2)) + 0.5 * alpha + (alpha + 1.0)
    }

    /// 遷移確率(中間)を返します。
    fn calc_prob_mid(alpha: f64, rate_fluctuation_var: f64, next_rate_fluctuation: f64) -> f64 {
        1.0 - rate_fluctuation_var / (2.0 * next_rate_fluctuation.powi(2)) - alpha.powi(2)
    }

    /// 遷移確率(下落)を返します。
    fn calc_prob_down(alpha: f64, rate_fluctuation_var: f64, next_rate_fluctuation: f64) -> f64 {
        rate_fluctuation_var / (2.0 * next_rate_fluctuation.powi(2)) + 0.5 * alpha + (alpha - 1.0)
    }

    /// 遷移先(上昇)のインデックスを返します。
    pub fn get_transition_index_up(&self) -> usize {
        self.transition_to[2]
    }

    /// 遷移先(中間)のインデックスを返します。
    pub fn get_transition_index_mid(&self) -> usize {
        self.transition_to[1]
    }

    /// 遷移先(下落)のインデックスを返します。
    pub fn get_transition_index_down(&self) -> usize {
        self.transition_to[0]
    }

    // Arrow Debreu price
}
