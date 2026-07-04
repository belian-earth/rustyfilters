//! Guided filter (He et al. 2013): the `rf_guided()` backend.
//!
//! A local linear model `out = a * I + b` is fitted per window between the
//! guide `I` and the input `p` (`a = cov(I, p) / (var(I) + eps)`,
//! `b = mean(p) - a * mean(I)`), then the coefficients are box-averaged and
//! applied at each cell. Everything reduces to box means on the O(1)
//! moments path, so the cost is independent of window size. Self-guided
//! (`I = p`) gives an edge-preserving smoother.

use crate::engine::{Dims, Edge, Win};
use crate::focal::{check_geom, fill_out, parse_edge, run_stat_m};
use extendr_api::prelude::*;

const NAN: f64 = f64::NAN;

fn box_mean(src: &[f64], d: Dims, w: Win, edge: Edge, na_omit: bool, out: &mut [f64]) {
    run_stat_m(
        src,
        d,
        w,
        edge,
        na_omit,
        |_, m| if m.n == 0 { NAN } else { m.mean() },
        out,
    );
}

#[allow(clippy::too_many_arguments)]
fn guided_impl(
    p: &[f64],
    guide: &[f64],
    d: Dims,
    w: Win,
    eps: f64,
    edge: Edge,
    na_omit: bool,
    out: &mut [f64],
) {
    let n = p.len();
    // Mask the union of missing cells into both inputs so every box mean
    // sees the same validity pattern (only matters for cross-guidance).
    let (i_m, p_m): (Vec<f64>, Vec<f64>) = guide
        .iter()
        .zip(p)
        .map(|(&i, &v)| {
            if i.is_nan() || v.is_nan() {
                (NAN, NAN)
            } else {
                (i, v)
            }
        })
        .unzip();
    let mut m_i = vec![0.0; n];
    let mut m_p = vec![0.0; n];
    let mut m_ii = vec![0.0; n];
    let mut m_ip = vec![0.0; n];
    box_mean(&i_m, d, w, edge, na_omit, &mut m_i);
    box_mean(&p_m, d, w, edge, na_omit, &mut m_p);
    let prod: Vec<f64> = i_m.iter().map(|v| v * v).collect();
    box_mean(&prod, d, w, edge, na_omit, &mut m_ii);
    let prod: Vec<f64> = i_m.iter().zip(&p_m).map(|(a, b)| a * b).collect();
    box_mean(&prod, d, w, edge, na_omit, &mut m_ip);
    drop(prod);
    // a and b per window, reusing the mean buffers.
    let mut a = m_ip;
    let mut b = m_p;
    for k in 0..n {
        let ak = (a[k] - m_i[k] * b[k]) / (m_ii[k] - m_i[k] * m_i[k] + eps);
        b[k] -= ak * m_i[k];
        a[k] = ak;
    }
    let (mut m_a, mut m_b) = (m_i, m_ii);
    box_mean(&a, d, w, edge, na_omit, &mut m_a);
    box_mean(&b, d, w, edge, na_omit, &mut m_b);
    for k in 0..n {
        out[k] = m_a[k] * guide[k] + m_b[k];
    }
}

/// Guided filter over a stack of column-major layers. `guide` must have the
/// same length as `x` (pass `x` itself for self-guidance).
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_guided_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    guide: &[f64],
    eps: f64,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    if guide.len() != x.len() {
        throw_r_error("`guide` must have the same dimensions as `x`");
    }
    if !eps.is_finite() || eps <= 0.0 {
        throw_r_error("`eps` must be a positive finite number");
    }
    let edge = parse_edge(edge, edge_value);
    fill_out(x.len(), |out| {
        guided_impl(x, guide, d, w, eps, edge, na_omit, out);
    })
}

extendr_module! {
    mod guided;
    fn rf_guided_rs;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_field_is_invariant() {
        let d = Dims { nr: 5, nc: 4 };
        let x = vec![3.0; 20];
        let mut out = vec![0.0; 20];
        guided_impl(&x, &x, d, Win { hr: 1, hc: 1 }, 1e-2, Edge::Shrink, true, &mut out);
        assert!(out.iter().all(|v| (v - 3.0).abs() < 1e-12));
    }

    #[test]
    fn large_eps_degenerates_towards_double_box_mean() {
        // As eps -> inf, a -> 0 and b -> mean(p), so out -> box_mean(mean(p)).
        let d = Dims { nr: 4, nc: 4 };
        let x: Vec<f64> = (0..16).map(|i| i as f64).collect();
        let w = Win { hr: 1, hc: 1 };
        let mut out = vec![0.0; 16];
        guided_impl(&x, &x, d, w, 1e12, Edge::Shrink, true, &mut out);
        let mut m1 = vec![0.0; 16];
        box_mean(&x, d, w, Edge::Shrink, true, &mut m1);
        let mut m2 = vec![0.0; 16];
        box_mean(&m1, d, w, Edge::Shrink, true, &mut m2);
        for (a, b) in out.iter().zip(&m2) {
            assert!((a - b).abs() < 1e-6);
        }
    }

    #[test]
    fn na_centre_stays_na() {
        let d = Dims { nr: 3, nc: 3 };
        let mut x = vec![1.0; 9];
        x[4] = NAN;
        let mut out = vec![0.0; 9];
        guided_impl(&x, &x, d, Win { hr: 1, hc: 1 }, 1e-2, Edge::Shrink, true, &mut out);
        assert!(out[4].is_nan());
        assert!(!out[0].is_nan());
    }
}
