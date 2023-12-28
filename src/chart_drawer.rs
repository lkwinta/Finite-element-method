use std::ops::Range;
use plotters::backend::BitMapBackend;
use plotters::chart::ChartBuilder;
use plotters::prelude::*;

use crate::function::Function;

pub fn draw_chart_single(title: &str, filename: &str, x_spec: Range<f32>, y_spec: Range<f32>, f: &Function) {
    let filename = "plots/".to_owned() + filename;
    let chart_root = BitMapBackend::new(&filename, (1920, 1080)).into_drawing_area();
    chart_root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&chart_root)
        .caption(title, ("times-new-roman", 20).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(x_spec, y_spec)
        .unwrap();
    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            (0..=200).map(|x| x as f32 * (0.01)).map(|x| (x, f.get_value(x as f64) as f32)),
            ShapeStyle::from(&RED).stroke_width(4),
        ))
        .unwrap();

    chart_root.present().unwrap();
}

pub fn draw_chart_array(title: &str, filename: &str, x_spec: Range<f32>, y_spec: Range<f32>, f_array: &Vec<Function>) {
    let filename = "plots/".to_owned() + filename;
    let chart_root = BitMapBackend::new(&filename, (1920, 1080)).into_drawing_area();
    chart_root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&chart_root)
        .caption(title, ("times-new-roman", 20).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(x_spec, y_spec)
        .unwrap();
    chart.configure_mesh().draw().unwrap();

    for (k, f) in f_array.iter().enumerate() {
        chart
            .draw_series(LineSeries::new(
                (0..=200).map(|x| x as f32 * (0.01)).map(|x| (x, f.get_value(x as f64) as f32)),
                ShapeStyle::from(
                    &RGBColor(
                        (k*130 % 255) as u8,
                        (k*190 % 255) as u8,
                        (k*230 % 255) as u8
                    )).stroke_width(4),
            ))
            .unwrap();
    }

    chart_root.present().unwrap();
}

