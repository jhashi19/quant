use super::curve::{self, Curve};
use super::hw_lib::HullWhite;
use super::node::Node;

#[derive(Clone, Debug)]
pub struct Tree {
    pub hw: HullWhite,              // Hull-White モデルのパラメータを持つインスタンス
    pub time_vec: Vec<f64>,         // 時間方向のグリッドのベクトル
    pub rate_interval: Vec<f64>,    // 金利方向のグリッドの間隔
    pub adjusting_params: Vec<f64>, // 各Nodeの金利をマーケットにadjustさせるパラメータ
    pub tree: Vec<Vec<Node>>,       // Treeの本体
}

// TODO Tree::new でTreeの構築が完了するようにロジックを修正する。
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

        let init_tree = Self {
            hw,
            time_vec,
            rate_interval,
            adjusting_params: vec![],
            tree: vec![vec![]],
        };

        let base_tree = init_tree.construct_base_tree();
        let adjusted_tree = base_tree.adjust_tree();
        adjusted_tree
    }

    /// branching process を構築します。
    fn construct_base_tree(self) -> Tree {
        //センタリングされたインデックス(-n～+n)と実際のVectorのインデックス(0～2n+1)に注意。
        let tree_length = self.time_vec.len() - 1; // time_vecの最後の要素の時刻のNodeは不要なため-1する。
        let mut adjusting_params: Vec<f64> = Vec::with_capacity(tree_length);
        let mut tree: Vec<Vec<Node>> = Vec::with_capacity(tree_length);
        tree[0] = Vec::with_capacity(1);
        // 最初のNodeは前の時刻が存在しないため、別で計算が必要。
        let a = self.get_a(self.time_vec[0]);
        let sigma = self.get_sigma(self.time_vec[0]);
        let rate_interval = self.rate_interval[0];
        let next_rate_interval = self.rate_interval[1]; // 次の時刻の金利方向の変動幅
        let time_interval = self.time_vec[1] - self.time_vec[0];
        let node = Self::create_node(
            a,
            sigma,
            0,
            0,
            rate_interval,
            next_rate_interval,
            time_interval,
            &Vec::new(),
            0.0,
            0.0,
        );
        tree[0][0] = node;
        let df = curve::df(Curve::Ois, self.time_vec[1]);
        adjusting_params[0] = Self::calc_adjusting_param(&tree[0], time_interval, df);

        // ２つめ以降のNode
        for i in 1..tree_length - 1 {
            let a = self.get_a(self.time_vec[i]);
            let sigma = self.get_sigma(self.time_vec[i]);

            // １つ前の時刻のNodeの金利方向のインデックスの最大値
            let previous_max_rate_index = tree[i - 1].len() - 1;
            // 今回のループで構築するNodeの金利方向のインデックスの最大値 = 前の時刻のNodeの金利方向の最大のインデックスの遷移(上昇)
            let max_rate_index = tree[i - 1][previous_max_rate_index].get_transition_index_up();
            let rate_interval = self.rate_interval[i];
            let next_rate_interval = self.rate_interval[i + 1]; // 次の時刻の金利方向の変動幅
            let time_interval = self.time_vec[i + 1] - self.time_vec[i];
            let node_num = 2 * max_rate_index + 1;
            tree[i] = Vec::with_capacity(node_num as usize); // node_numは負になることはない。

            for j in 0..2 * max_rate_index + 1 {
                let p_node = &tree[i - 1];
                let p_time_interval = self.time_vec[i] - self.time_vec[i - 1];
                let adjusting_param = adjusting_params[i - 1];
                let node = Self::create_node(
                    a,
                    sigma,
                    j,
                    max_rate_index,
                    rate_interval,
                    next_rate_interval,
                    time_interval,
                    p_node,
                    p_time_interval,
                    adjusting_param,
                );
                tree[i][j as usize] = node;
            }
            let df = curve::df(Curve::Ois, self.time_vec[i + 1]);
            adjusting_params[i] = Self::calc_adjusting_param(&tree[i], time_interval, df);
        }

        Tree {
            hw: self.hw,
            time_vec: self.time_vec,
            rate_interval: self.rate_interval,
            adjusting_params,
            tree,
        }
    }

    /// Treeをマーケットデータにadjustさせます。
    /// 各Nodeにその時刻のadjusting_paramsを加算してマーケットに一致させます。
    fn adjust_tree(self) -> Tree {
        let mut tree: Tree = self;
        for i in 0..tree.tree.len() {
            for j in 0..tree.tree[i].len() {
                tree.tree[i][j].rate += tree.adjusting_params[i];
            }
        }
        tree
    }

    /// 指定したポイントでのpiecewise-constantなaを取得します。
    /// * `target` - 戻り値のポイント
    fn get_a(&self, target: f64) -> f64 {
        Self::get_piecewise_constant_value(&self.hw.a, &self.hw.a_interval, target)
    }

    /// 指定したポイントでのpiecewise-constantなsigmaを取得します。
    /// * `target` - 戻り値のポイント
    fn get_sigma(&self, target: f64) -> f64 {
        Self::get_piecewise_constant_value(&self.hw.sigma, &self.hw.sigma_interval, target)
    }

    /// piecewise-constantなパラメータの指定したポイントでの値を返します。<br>
    /// intervalの外側の値については端の値に一致するものとします。
    /// * `val` - piecewise-constantなパラメータの値のベクタ
    /// * `interval` - piecewise-constantなパラメータの値に対応する時間間隔のベクタ
    /// * `target` - 戻り値のポイント
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
    /// * `sigma` - Hull Whiteモデルパラメータのσ
    /// * `time_interval` - 次の時刻のNodeまでの時間間隔
    fn calc_rate_interval(sigma: f64, time_interval: f64) -> f64 {
        sigma * (3.0 * time_interval).powf(0.5)
    }

    /// Nodeの各パラメータを生成してそれらの値をもとにNodeのインスタンスを生成します。
    /// * `a` - Hull Whiteモデルパラメータのa
    /// * `sigma` - Hull Whiteモデルパラメータのσ
    /// * `index` - 同一時刻のNodeの金利方向の実際のインデックス
    /// * `max_rate_index` - Nodeの金利方向のインデックスの最大値
    /// * `rate_interval` - Nodeの金利方向の間隔
    /// * `next_rate_interval` - 次の時刻のNodeの金利方向の間隔
    /// * `time_interval` - 次の時刻までの時間間隔
    /// * `p_nodes` - 前の時刻のNodeのベクタ
    /// * `p_time_interval` - 前の時刻からの時間間隔
    /// * `p_adjusting_param` - 前の時刻の金利の調整項
    fn create_node(
        a: f64,
        sigma: f64,
        index: isize,
        max_rate_index: isize,
        rate_interval: f64,
        next_rate_interval: f64,
        time_interval: f64,
        p_nodes: &Vec<Node>,
        p_time_interval: f64,
        p_adjusting_param: f64,
    ) -> Node {
        // センタリングされたインデックス
        let centering_index = index - max_rate_index;
        let rate = centering_index as f64 * rate_interval;
        let rate_fluctuation_mean = Node::calc_fluctuation_mean(a, rate, time_interval);
        let rate_fluctuation_var = Node::calc_fluctuation_var(a, sigma, time_interval);

        // 遷移先Node(中間)のセンタリングされたインデックス
        let transition_to_mid_index = ((rate + rate_fluctuation_mean) / next_rate_interval).round();
        let next_node_mid_rate = transition_to_mid_index * next_rate_interval;
        let (prob_up, prob_mid, prob_down) = Node::calc_transition_prob(
            rate,
            rate_fluctuation_mean,
            rate_fluctuation_var,
            next_node_mid_rate,
            next_rate_interval,
        );

        let arrow_debreu =
            Self::calc_arrow_debreu(p_nodes, p_time_interval, p_adjusting_param, centering_index);

        Node::new(
            rate,
            arrow_debreu,
            transition_to_mid_index as isize,
            rate_fluctuation_mean,
            rate_fluctuation_var,
            prob_up,
            prob_mid,
            prob_down,
        )
    }

    /// adjusting_paramsの各値を返します。（マーケットの金利とベースのツリーの金利を一致させるための調整項）
    /// * `nodes` - Nodeのベクタ
    /// * `time_interval` - 次の時刻までの時間間隔
    /// * `df` - マーケットのレートから計算したディスカウントファクター
    fn calc_adjusting_param(nodes: &Vec<Node>, time_interval: f64, df: f64) -> f64 {
        let mut unadjusted_df = 0.0;
        for node in nodes {
            unadjusted_df += node.arrow_debreu * (-node.rate * time_interval).exp()
        }
        (unadjusted_df.ln() - df.ln()) / time_interval
    }

    /// 各Nodeに設定するarrow_debreuの値を返します。
    /// 計算に使用する各パラメータは1つ前の時刻のものを使用するためp_から始まる（previous）
    /// * `p_nodes` - 前の時刻のNodeのベクタ
    /// * `p_time_interval` - 前の時刻からの時間間隔
    /// * `p_adjusting_param` - 前の時刻の金利の調整項
    /// * `centering_index` - arrow_debreu を計算するNodeのセンタリングされたインデックス
    fn calc_arrow_debreu(
        p_nodes: &Vec<Node>,
        p_time_interval: f64,
        p_adjusting_param: f64,
        centering_index: isize,
    ) -> f64 {
        // TODO 最初のNode用の分岐。もっといい方法がないか検討する。
        if p_nodes.len() == 0 {
            return 1.0;
        }
        let mut arrow_debreu = 0.0;
        for p_node in p_nodes {
            // TODO ループの中で以下の条件に掛からないものが多いので効率的にforを抜ける方法を検討する
            let calc_arrow_debreu_component = |transition_prob: f64| -> f64 {
                p_node.arrow_debreu
                    * transition_prob
                    * (-(p_adjusting_param + p_node.rate) * p_time_interval).exp()
            };
            if p_node.get_transition_index_up() == centering_index {
                arrow_debreu += calc_arrow_debreu_component(p_node.prob_up);
            }
            if p_node.get_transition_index_mid() == centering_index {
                arrow_debreu += calc_arrow_debreu_component(p_node.prob_mid);
            }
            if p_node.get_transition_index_down() == centering_index {
                arrow_debreu += calc_arrow_debreu_component(p_node.prob_down);
            }
        }
        arrow_debreu
    }

    // 各種商品をバックワードでプライシング → ファイル分ける
}
