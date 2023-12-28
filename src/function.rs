use std::slice::Iter;

pub struct Function {
    functions: Vec<((f64, f64), Box<dyn Fn(f64) -> f64>)>,
    other_case_function: Box<dyn Fn(f64) -> f64>
}

impl Function {
    pub fn new<F: Fn(f64) -> f64 + 'static>(base_case: F) -> Self {
        Self {
            functions: Vec::new(),
            other_case_function: Box::new(base_case)
        }
    }

    pub fn add_case<F: Fn(f64) -> f64 + 'static>(mut self, range: (f64, f64), function: F) -> Self{
        self.functions.push((range, Box::new(function)));
        self
    }

    pub fn get_value(&self, x: f64) -> f64{
        for (range, function) in &self.functions {
            if range.0 <= x && range.1 > x {
                return function(x);
            }
        }
        (self.other_case_function)(x)
    }

    pub fn cases_iter(&self) -> Iter<'_, ((f64, f64), Box<dyn Fn(f64) -> f64>)> {
        self.functions.iter()
    }
}