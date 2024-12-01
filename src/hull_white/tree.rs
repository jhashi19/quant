use super::node::Node;

#[derive(Clone, Debug)]
pub struct Tree {
    pub time_vec: Vec<f64>,      // 時間方向のグリッドのベクトル
    pub rate_interval: Vec<f64>, // 金利方向のグリッドの間隔
    pub tree: Vec<Vec<Node>>,    // Treeの本体
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
