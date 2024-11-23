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

#[derive(Clone, Debug)]
pub struct Tree {
    pub time_vec: Vec<f64>,      // 時間方向のグリッドのベクトル
    pub rate_interval: Vec<f64>, // 金利方向のグリッドの間隔
    pub tree: Vec<Vec<Node>>,    // Treeの本体
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

impl Tree {
    // 初期処理
    // 時間方向のベクトルと金利方向の間隔を設定する。
    pub fn new(sigma: &f64, time_vec: Vec<f64>) -> Self {
        let interval_num = time_vec.len() - 1;
        let mut rate_interval = Vec::with_capacity(interval_num);
        for i in 0..interval_num {
            let time_interval = time_vec[i + 1] - time_vec[i];
            rate_interval[i] = Self::calc_rate_interval(sigma, time_interval);
        }
        Self {
            time_vec,
            rate_interval,
            tree: vec![vec![]],
        }
    }

    /// branching process を構築する。
    // pub fn construct_base_tree(self) -> Tree {}

    /// Treeをマーケットデータにadjustさせる。
    // pub fn adjust_tree(self) -> {}

    /// 金利方向のグリッドの間隔を返します。
    fn calc_rate_interval(sigma: &f64, time_interval: f64) -> f64 {
        sigma * (3.0 * time_interval).powf(0.5)
    }

    // Arrow-Debreu Tree作る？

    // 各種商品をバックワードでプライシング → ファイル分ける
}
