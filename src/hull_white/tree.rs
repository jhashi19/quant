pub struct Tree {
    pub delta_t: f64,
    pub delta_x: f64,
    pub mean: f64,
    pub var: f64,
    pub x_idx_max: usize,
}

impl Tree {
    // 初期処理
    // input: a, sigma, delta_t(ノードの時間間隔(横軸)), h(金利の間隔(縦軸)を計算するときのinput)
    // output: 金利の間隔の期待値, 金利の間隔の分散, 金利の間隔(縦軸)
    pub fn new(a: f64, sigma: f64, delta_t: f64, h: f64) -> Self {
        let mean = (-a * delta_t).exp() - 1.0;
        let var = sigma.powi(2) * 0.5 / a * (1.0 - (-2.0 * a * delta_t).exp());
        let x_idx_max = ((1.0 - (1.0 - 1.0 / h).sqrt()) / mean).abs().ceil() as usize * 2;
        let delta_x = (h * var).sqrt();
        Tree {
            delta_t,
            delta_x,
            mean,
            var,
            x_idx_max,
        }
    }

    // theta=0 とした場合のツリーの構築

    // 上記のツリーのイールドカーブへのフィッティング

    // 各種商品をバックワードでプライシング
}
