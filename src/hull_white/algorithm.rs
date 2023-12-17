pub struct Newton {
    func: fn(x: f64) -> f64, // ターゲットとなる1次元の関数。 =0の形であらわした左辺。
    func_deriv: fn(x: f64) -> f64, // funcの導関数
    x_max: f64,              // funcの表す方程式の解が取りうる値の範囲の上限
    x_min: f64,              // funcの表す方程式の解が取りうる値の範囲の下限
    threshold: f64,          // 収束判定の閾値
}

impl Newton {
    pub fn new(
        func: fn(x: f64) -> f64,
        func_deriv: fn(x: f64) -> f64,
        x_max: f64,
        x_min: f64,
        threshold: f64,
    ) -> Self {
        Newton {
            func,
            func_deriv,
            x_max,
            x_min,
            threshold,
        }
    }

    pub fn find_root(&self) -> f64 {
        let max_iter = 20;
        let mut root = 0.5 * (self.x_max + self.x_min);
        for _ in 0..max_iter {
            let f = (self.func)(root);
            let df = (self.func_deriv)(root);
            let dx = f / df;
            root -= dx;
            if (self.x_max - root) * (root - self.x_min) < 0.0 {
                panic!("Root is out of range [{},{}]", self.x_min, self.x_max);
            }
            if dx.abs() < self.threshold {
                return root;
            }
        }
        panic!("Maximun number of iterations exceeded");
    }

    pub fn find_root_safe(&self) -> f64 {
        let max_iter = 100;
        let mut xf_high: f64; // x_maxとx_minのうち、funcが大きい値をとる方
        let mut xf_low: f64; // x_maxとx_minのうち、funcが小さい値をとる方
        let fr = (self.func)(self.x_max);
        let fl = (self.func)(self.x_min);
        if (fl > 0.0 && fr > 0.0) || (fl < 0.0 && fr < 0.0) {
            panic!("Root must be bracketed in [{},{}]", self.x_min, self.x_max);
        }
        if fl == 0.0 {
            return self.x_min;
        }
        if fr == 0.0 {
            return self.x_max;
        }
        if fl < 0.0 {
            xf_high = self.x_max;
            xf_low = self.x_min;
        } else {
            xf_high = self.x_min;
            xf_low = self.x_max;
        }
        let mut root = 0.5 * (self.x_max + self.x_min);
        let mut dxold = self.x_max - self.x_min;
        let mut dx = dxold;
        let mut f = (self.func)(root);
        let mut df = (self.func_deriv)(root);
        for _ in 0..max_iter {
            if ((root - xf_high) * df - f) * ((root - xf_low) * df - f) > 0.0
                || (2.0 * f).abs() > (dxold * df).abs()
            {
                dxold = dx;
                dx = 0.5 * (xf_high - xf_low);
                root = xf_low + dx;
            } else {
                dxold = dx;
                dx = f / df;
                root -= dx;
            }
            if dx.abs() < self.threshold {
                return root;
            }
            f = (self.func)(root);
            df = (self.func_deriv)(root);
            if f < 0.0 {
                xf_low = root;
            } else {
                xf_high = root;
            }
        }
        panic!("Maximun number of iterations exceeded");
    }
}

pub struct LevenbergMarquardt {}

impl LevenbergMarquardt {
    pub fn new() -> Self {
        LevenbergMarquardt {}
    }

    pub fn fit() {}
}

#[cfg(test)]
mod tests {
    use super::*;

    fn func(x: f64) -> f64 {
        x.powi(2) - 2.0
    }
    fn func_deriv(x: f64) -> f64 {
        2.0 * x
    }

    #[test]
    fn test_find_root() {
        let ins = Newton::new(func, func_deriv, 2.0, 0.0, 1.0e-7);
        let root = ins.find_root();
        assert!((root - 2_f64.powf(0.5)).abs() < 1.0e-10);
    }

    #[test]
    fn test_find_root_safe() {
        let ins = Newton::new(func, func_deriv, 2.0, 0.0, 1.0e-7);
        let root = ins.find_root_safe();
        assert!((root - 2_f64.powf(0.5)).abs() < 1.0e-10);
    }
}
