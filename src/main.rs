#![allow(uncommon_codepoints)]
#![allow(mixed_script_confusables)]
#![allow(non_snake_case)]

mod chart_drawer;
mod function;

use gauss_quad::GaussLegendre;
use nalgebra::{DMatrix, DVector};

use Fn;
use function::Function;

fn get_e_k(x_k: f64, ε: f64) -> Function {
    let lower = x_k - ε;
    let upper = x_k + ε;
    Function::new(|_| 0.0)
        .add_case((lower, x_k), move |x| (x - lower) / (x_k - lower))
        .add_case((x_k, upper), move |x| (upper - x) / (upper - x_k))
}

fn get_de_k(x_k: f64, ε: f64) -> Function {
    let lower = x_k - ε;
    let upper = x_k + ε;
    Function::new(|_| 0.0)
        .add_case((lower, x_k), move |_| 1.0 / (x_k - lower))
        .add_case((x_k, upper), move |_| -1.0 / (upper - x_k))
}

fn integrate<F: Fn(f64) -> f64>(lower_bound: f64, upper_bound: f64, f: F) -> f64{
    let quad = GaussLegendre::init(3);
    quad.integrate(lower_bound, upper_bound, f)
}

fn B(u: &Function, du: &Function, ϕ: &Function, dϕ: &Function) -> f64 {
    let mut result = -u.get_value(0.0)*ϕ.get_value(0.0);

    for ((du_start, du_end), f_du) in du.cases_iter() {
        for ((dϕ_start, dϕ_end), f_dϕ) in dϕ.cases_iter() {
            if du_start > dϕ_end || dϕ_start > du_end {
                continue
            }
            let (int_start, int_end) =
                if dϕ_start >= du_start && dϕ_end <= du_end { (dϕ_start, dϕ_end) }
                else if dϕ_start < du_start { (du_start, dϕ_end) }
                else { (dϕ_start, du_end) };

            result +=
               if *int_end <= 1.0 { integrate(f64::max(0.0, *int_start), *int_end, |x| f_du(x) * f_dϕ(x)) }
               else if *int_start >= 1.0 { integrate(*int_start, f64::min(2.0, *int_end), |x| 2.0 * x * f_du(x) * f_dϕ(x)) }
               else {
                   integrate(f64::max(0.0, *int_start), 1.0, |x| f_du(x) * f_dϕ(x)) +
                   integrate(1.0,f64::min(2.0, *int_end), |x| 2.0 * x * f_du(x) * f_dϕ(x))
               }
        }
    }

    result
}

fn L(ϕ: &Function) -> f64 {
    let mut result = -20.0 * ϕ.get_value(0.0);

    for ((start, end), f_ϕ) in ϕ.cases_iter(){
        result += integrate(f64::max(0.0, *start), f64::min(2.0, *end), |x| 100.0*x*f_ϕ(x))
    }

    result
}


fn main() {
    const N: usize = 20;
    let a = 0.0;
    let b = 2.0;
    let Δ = (b - a)/N as f64;

    let mut B_matrix: DMatrix<f64> = nalgebra::DMatrix::zeros(N, N);
    let mut L_vector: DVector<f64> = nalgebra::DVector::zeros(N);
    let test_functions = (0..N).map(|k| get_e_k((k as f64)*Δ, Δ)).collect::<Vec<_>>();
    let derivatives_test_functions = (0..N).map(|k| get_de_k((k as f64)*Δ, Δ)).collect::<Vec<_>>();

    for j in 0..N {
        for i in 0..N {
            let element = B_matrix.get_mut((j, i)).unwrap();
            *element = B( &test_functions[i], &derivatives_test_functions[i], &test_functions[j], &derivatives_test_functions[j]);
        }
        L_vector[j] = L( &test_functions[j]);
    }

    let W_vector = B_matrix.lu().solve(&L_vector).unwrap();

    chart_drawer::draw_chart_array("Funkcje testujące",
                                   "test_functions.png",
                                   -0.1f32..2.1f32, -0.1f32..1.5f32,
                                   &test_functions);

    let u =
        Function::new(move |x| W_vector.iter().zip(&test_functions).map(|(w_k, e_k)| w_k*e_k.get_value(x)).sum::<f64>() );

    println!("u(0) = {}", u.get_value(0.0));
    let h = 10e-12;
    println!("du/dx(0) = {}", (u.get_value(0.0 + h) - u.get_value(0.0))/h);
    println!("u(2) = {}", u.get_value(2.0));

    chart_drawer::draw_chart_single("Wynik rozwiązania równania",
                                    "equation_result.png",
                                    -0.1f32..2.1f32, -80.0f32..15.0f32,
                                    &u);

}
