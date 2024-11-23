// ■導出の流れ
// ①連立方程式を導出する。
// ②Thomas法で連立方程式を解く。
// ③連立方程式の解から3次多項式の係数を算出する。
// ④targetにおけるレートを算出する。

pub fn cubic_spline(dates: &Vec<f64>, rates: &Vec<f64>, target: &f64) -> f64 {
    // TODO 引数チェックを入れる

    let equation_matrix = set_up_equation(&dates, &rates);
    let mut solution_of_equation = solve_equation(equation_matrix);

    // natural spline
    let mut second_derive = vec![0.0];
    second_derive.append(&mut solution_of_equation);
    second_derive.append(&mut vec![0.0]);

    // targetの値が属する区間の始点側のindexを算出する。
    // extrapolationとなる点は分岐に入る頻度が少ないので後ろに記述
    let mut target_idx = 0;
    if &dates[0] <= target && target <= &dates[dates.len() - 1] {
        for (i, date) in dates.iter().enumerate() {
            if date > &target {
                target_idx = i - 1;
                break;
            }
        }
    } else if target <= &dates[0] {
        target_idx = 0;
    } else {
        target_idx = dates.len() - 2;
    }

    // 多項式の係数を算出する。a_j * (x - x_j)^3 + b_j * (x - x_j)^2 + c_j * (x - x_j) + d_j
    let a = (second_derive[target_idx + 1] - second_derive[target_idx])
        / (6.0 * (dates[target_idx + 1] - dates[target_idx]));
    let b = second_derive[target_idx] * 0.5;
    let c = (rates[target_idx + 1] - rates[target_idx])
        / (dates[target_idx + 1] - dates[target_idx]) as f64
        - (dates[target_idx + 1] - dates[target_idx]) as f64
            * (2.0 * second_derive[target_idx] + second_derive[target_idx + 1])
            / 6.0;
    let d = rates[target_idx];

    let period = (target - dates[target_idx]) as f64;
    a * period.powf(3.0) + b * period.powf(2.0) + c * period + d
}

// cubic spline のための連立方程式を表す3行対角行列の上中下の対角成分で構成される各ベクタと右辺のベクタをtupleで返す
fn set_up_equation(dates: &Vec<f64>, rates: &Vec<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let dates_num = dates.len();

    let mut period: Vec<f64> = Vec::with_capacity(dates_num - 1);
    for i in 0..dates_num - 1 {
        period.push(dates[i + 1] - dates[i]);
    }

    let mut upper_diagonal: Vec<f64> = Vec::with_capacity(dates_num - 2);
    let mut diagonal: Vec<f64> = Vec::with_capacity(dates_num - 2);
    let mut lower_diagonal: Vec<f64> = Vec::with_capacity(dates_num - 2);
    let mut rs_vector: Vec<f64> = Vec::with_capacity(dates_num - 2);

    for i in 0..dates_num - 2 {
        upper_diagonal.push(period[i + 1] as f64);
        diagonal.push(2.0 * (dates[i + 2] - dates[i]) as f64);
        lower_diagonal.push(period[i] as f64);
        rs_vector.push(
            6.0 * ((rates[i + 2] - rates[i + 1]) / period[i + 1] as f64
                - (rates[i + 1] - rates[i]) / period[i] as f64),
        );
    }
    (upper_diagonal, diagonal, lower_diagonal, rs_vector)
}

// Thomas algorithm(Tridiagonal matrix algorithm)で3行対角行列で表される連立方程式の解を求めて返す。
fn solve_equation(
    (mut upper_diagonal, mut diagonal, lower_diagonal, mut rs_vector): (
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
        Vec<f64>,
    ),
) -> Vec<f64> {
    upper_diagonal[0] /= diagonal[0];
    rs_vector[0] /= diagonal[0];
    diagonal[0] = 1.0;
    for i in 1..upper_diagonal.len() {
        upper_diagonal[i] /= diagonal[i] - lower_diagonal[i] * upper_diagonal[i - 1];
        rs_vector[i] = (rs_vector[i] - rs_vector[i - 1] * lower_diagonal[i])
            / (diagonal[i] - upper_diagonal[i - 1] * lower_diagonal[i]);
        diagonal[i] = 1.0;
    }

    for i in (0..rs_vector.len() - 1).rev() {
        rs_vector[i] = rs_vector[i] - upper_diagonal[i] * rs_vector[i + 1];
    }
    rs_vector
}

#[cfg(test)]
mod tests {
    use crate::hull_white::interpolation::cubic_spline;

    #[test]
    fn test_cubic_spline() {
        let dates = vec![
            100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
        ];
        let rates = vec![
            0.889, 1.085664, 1.05, 1.12, 1.327, 2.0512, 2.578, 2.2245, 2.115143, 2.356,
        ];
        let result = cubic_spline(&dates, &rates, &118.0);
        assert_eq!(result, 0.936556517871382);
    }
}
