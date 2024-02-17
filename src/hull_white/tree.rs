pub struct Tree {
    pub delta_t: f64,          // Treeの時間の区間幅
    pub delta_x: f64,          // Treeの金利の区間幅
    pub mean: f64,             // 金利のdelta_tの間における変化の平均
    pub var: f64,              // 金利のdelta_tの間における変化の分散
    pub x_idx_max: usize,      // Treeの金利方向の最大ノードのインデックス
    pub p_norm: Vec<Vec<f64>>, // 金利の通常の推移確率
    pub p_max: Vec<f64>,       // 金利方向の最大ノードの推移確率
    pub p_min: Vec<f64>,       // 金利方向の最小ノードの推移確率
    pub base: Vec<Vec<f64>>,   // theta=0としたときのTree
}

impl Tree {
    // 初期処理
    // input: a, sigma, delta_t(ノードの時間間隔(横軸)), h(金利の間隔(縦軸)を計算するときのinput)
    // output: 金利の間隔の期待値, 金利の間隔の分散, 金利の間隔(縦軸)
    // TODO インプットとアウトプットとして何が必要か整理できていないため整理する。
    // TODO サンプルの実装なので、期間は10年固定でツリーを実装する、とか決めてしまってもよいかも。細かい実装の柔軟性は本質的ではない。
    pub fn new(a: f64, sigma: f64, delta_t: f64, h: f64) -> Self {
        let mean = (-a * delta_t).exp() - 1.0;
        let var = sigma.powi(2) * 0.5 / a * (1.0 - (-2.0 * a * delta_t).exp());
        let x_idx_max = ((1.0 - (1.0 - 1.0 / h).sqrt()) / mean).abs().ceil() as usize * 2;
        // p_normは(i, 0)がup、(i, 1)がmiddle、(i, 2)がdown
        // インデックスを -n～n → 0～2n+1 に変換する必要がある。
        let mut p_norm = Vec::with_capacity(2 * x_idx_max + 1);
        p_norm[x_idx_max + 1] = vec![0.5 * (1.0 / h), 1.0 - 1.0 / h, 0.5 * (1.0 / h)];
        for i in 1..=x_idx_max {
            // プラス側
            let idx_mean = i as f64 * mean;
            let prob_vec = vec![
                0.5 * (1.0 / h + idx_mean.powi(2) + idx_mean),
                1.0 - 1.0 / h - idx_mean.powi(2),
                0.5 * (1.0 / h + idx_mean.powi(2) - idx_mean),
            ];
            p_norm[x_idx_max + 1 + i] = prob_vec;
            // マイナス側
            let idx_mean_minus = -1.0 * i as f64 * mean;
            let prob_vec = vec![
                0.5 * (1.0 / h + idx_mean_minus.powi(2) + idx_mean_minus),
                1.0 - 1.0 / h - idx_mean_minus.powi(2),
                0.5 * (1.0 / h + idx_mean_minus.powi(2) - idx_mean_minus),
            ];
            p_norm[x_idx_max + 1 - i] = prob_vec;
        }

        let idx_max_mean = x_idx_max as f64 * mean;
        // p_maxは(up, middle, down)
        let p_max = vec![
            1.0 + 0.5 * (1.0 / h + idx_max_mean.powi(2) + 3.0 * idx_max_mean),
            -(1.0 / h + idx_max_mean.powi(2) + 2.0 * idx_max_mean),
            0.5 * (1.0 / h + idx_max_mean.powi(2) + idx_max_mean),
        ];
        // p_minは(up, middle, down)
        let p_min = vec![
            0.5 * (1.0 / h + idx_max_mean.powi(2) + idx_max_mean),
            -(1.0 / h + idx_max_mean.powi(2) + 2.0 * idx_max_mean),
            1.0 + 0.5 * (1.0 / h + idx_max_mean.powi(2) + 3.0 * idx_max_mean),
        ];
        let delta_x = (h * var).sqrt();
        let base = Tree::construct_base_tree();
        Tree {
            delta_t,
            delta_x,
            mean,
            var,
            x_idx_max,
            p_norm,
            p_max,
            p_min,
            base,
        }
    }

    // サンプル
    fn construct_base_tree() -> Vec<Vec<f64>> {
        let mut base: Vec<Vec<f64>> = vec![vec![0.0; 10]; 10];
        base
    }

    // theta=0 とした場合のツリーの構築
    //   ノードごとの推移確率（Treeに持たせて初期処理で設定する方がよいかも）
    //   ↑そこまでやってしまえばtheta=0の場合のツリーの構築までを初期処理でやる方がよいかも

    // 上記のツリーのイールドカーブへのフィッティング

    // 各種商品をバックワードでプライシング
}
