use super::hw_lib::HullWhite;
use super::node::Node;

#[derive(Clone, Debug)]
pub struct Tree {
    pub hw: HullWhite,           // Hull-White モデルのパラメータを持つインスタンス
    pub time_vec: Vec<f64>,      // 時間方向のグリッドのベクトル
    pub rate_interval: Vec<f64>, // 金利方向のグリッドの間隔
    pub tree: Vec<Vec<Node>>,    // Treeの本体
}

impl Tree {
    // 初期処理
    // 時間方向のベクトルと金利方向の間隔を設定する。
    pub fn new(hw: HullWhite, time_vec: Vec<f64>) -> Self {
        let time_interval_num = time_vec.len();
        let mut rate_interval = Vec::with_capacity(time_interval_num);
        rate_interval[0] = 0.0; // 0番目はTreeのスタートでNodeが1個なのでrate_intervalは使用しないため0としておく。
        for i in 1..time_interval_num {
            let time_interval = time_vec[i] - time_vec[i - 1];
            let sigma =
                Self::get_piecewise_constant_value(&hw.sigma, &hw.sigma_interval, time_vec[i]);
            rate_interval[i] = Self::calc_rate_interval(sigma, time_interval);
        }
        Self {
            hw,
            time_vec,
            rate_interval,
            tree: vec![vec![]],
        }
    }

    /// branching process を構築します。
    pub fn construct_base_tree(&self) -> Tree {
        //センタリングされたインデックス(-n～+n)と実際のVectorのインデックス(0～2n+1)に注意。
        let tree_length = self.time_vec.len();
        let mut tree: Vec<Vec<Node>> = Vec::with_capacity(tree_length);
        // 最初のNode とりあえず普通にNodeを構成する。
        let a = self.get_a(self.time_vec[0]);
        let sigma = self.get_sigma(self.time_vec[0]);
        let mut node = Node::new();
        // Nodeの各パラメータを計算する
        let time_interval = self.time_vec[1] - self.time_vec[0];
        let rate_fluctuation_mean = Node::calc_fluctuation_mean(&a, node.rate, time_interval);
        let rate_fluctuation_var = Node::calc_fluctuation_var(&a, &sigma, time_interval);
        // 遷移確率はほかのNodeの情報を使う
        let next_rate_fluctuation = self.rate_interval[1]; // 次の時刻の金利方向の変動幅
        let transition_to_mid_index =
            ((node.rate + rate_fluctuation_mean) / next_rate_fluctuation).round(); // 遷移先Node(中間)のセンタリングされたインデックス
        let next_node_mid_rate = transition_to_mid_index * next_rate_fluctuation;
        let (prob_up, prob_mid, prob_down) = Node::calc_transition_prob(
            node.rate,
            rate_fluctuation_mean,
            rate_fluctuation_var,
            next_node_mid_rate,
            next_rate_fluctuation,
        );
        let transition_to = vec![
            transition_to_mid_index as usize - 1,
            transition_to_mid_index as usize,
            transition_to_mid_index as usize + 1,
        ];

        node.transition_to = transition_to;
        // rateはインスタンス化した時の0.0のままでよい
        node.rate_fluctuation_mean = rate_fluctuation_mean;
        node.rate_fluctuation_var = rate_fluctuation_var;
        node.prob_up = prob_up;
        node.prob_mid = prob_mid;
        node.prob_down = prob_down;
        // arrow_debreuはadjust_treeで設定する。
        tree[0].push(node);

        // ２つめ以降のNode
        for i in 1..tree_length - 1 {
            let a = self.get_a(self.time_vec[i]);
            let sigma = self.get_sigma(self.time_vec[i]);

            // １つ前の時刻のNodeの金利方向のインデックスの最大値
            let previous_max_rate_index = tree[i - 1].len() - 1;
            // 今回のループで構築するNodeの金利方向のインデックスの最大値
            let max_rate_index = tree[i - 1][previous_max_rate_index].get_transition_index_up();
            for j in 0..2 * max_rate_index + 1 {
                // センタリングしたインデックスj
                let centering_j = j - max_rate_index;
                let mut node = Node::new();
                let rate = centering_j as f64 * self.rate_interval[i];
                let time_interval = self.time_vec[i] - self.time_vec[i - 1];
                let rate_fluctuation_mean = Node::calc_fluctuation_mean(&a, rate, time_interval);
                let rate_fluctuation_var = Node::calc_fluctuation_var(&a, &sigma, time_interval);

                // 遷移確率はほかのNodeの情報を使う
                let next_rate_fluctuation = self.rate_interval[i + 1]; // 次の時刻の金利方向の変動幅
                let transition_to_mid_index =
                    ((rate + rate_fluctuation_mean) / next_rate_fluctuation).round(); // 遷移先Node(中間)のセンタリングされたインデックス
                let next_node_mid_rate = transition_to_mid_index * next_rate_fluctuation;
                let (prob_up, prob_mid, prob_down) = Node::calc_transition_prob(
                    rate,
                    rate_fluctuation_mean,
                    rate_fluctuation_var,
                    next_node_mid_rate,
                    next_rate_fluctuation,
                );
                let transition_to = vec![
                    transition_to_mid_index as usize - 1,
                    transition_to_mid_index as usize,
                    transition_to_mid_index as usize + 1,
                ];

                node.transition_to = transition_to;
                node.rate = rate;
                node.rate_fluctuation_mean = rate_fluctuation_mean;
                node.rate_fluctuation_var = rate_fluctuation_var;
                node.prob_up = prob_up;
                node.prob_mid = prob_mid;
                node.prob_down = prob_down;
                // arrow_debreuはadjust_treeで設定する。
                tree[i].push(node);
            }
        }

        // 最後の時刻のNode
        // ここではrateだけでよい？
        let last_node_index = tree_length - 1;
        // １つ前の時刻のNodeの金利方向のインデックスの最大値
        let previous_max_rate_index = tree[last_node_index - 1].len() - 1;
        // 今回のループで構築するNodeの金利方向のインデックスの最大値
        let max_rate_index =
            tree[last_node_index - 1][previous_max_rate_index].get_transition_index_up();
        for j in 0..2 * max_rate_index + 1 {
            // センタリングしたインデックスj
            let centering_j = j - max_rate_index;
            let mut node = Node::new();
            let rate = centering_j as f64 * self.rate_interval[last_node_index];
            node.rate = rate;
            tree[last_node_index].push(node);
        }

        // cloneしなくてもよい方法はないか
        Tree {
            hw: self.hw.clone(),
            time_vec: self.time_vec.clone(),
            rate_interval: self.rate_interval.clone(),
            tree,
        }
    }

    /// Treeをマーケットデータにadjustさせます。
    // pub fn adjust_tree(self) -> Tree {}

    /// 指定したポイントでのpiecewise-constantなaを取得します。
    fn get_a(&self, target: f64) -> f64 {
        Self::get_piecewise_constant_value(&self.hw.a, &self.hw.a_interval, target)
    }

    /// 指定したポイントでのpiecewise-constantなsigmaを取得します。
    fn get_sigma(&self, target: f64) -> f64 {
        Self::get_piecewise_constant_value(&self.hw.sigma, &self.hw.sigma_interval, target)
    }

    /// piecewise-constantなパラメータの指定したポイントでの値を返します。<br>
    /// intervalの外側の値については端の値に一致するものとします。
    fn get_piecewise_constant_value(val: &Vec<f64>, interval: &Vec<f64>, target: f64) -> f64 {
        if target < interval[0] {
            return val[0];
        }
        if target > interval[interval.len() - 1] {
            return val[interval.len() - 1];
        }
        let mut val_index = 0;
        for i in 1..interval.len() {
            if target < interval[i] {
                break;
            }
            val_index += 1;
        }
        val[val_index]
    }

    /// 金利方向のグリッドの間隔を返します。
    fn calc_rate_interval(sigma: f64, time_interval: f64) -> f64 {
        sigma * (3.0 * time_interval).powf(0.5)
    }

    // Arrow-Debreu Tree作る？

    // 各種商品をバックワードでプライシング → ファイル分ける
}
