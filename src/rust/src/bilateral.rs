//! Bilateral filter: the `rf_bilateral()` backend.
//!
//! Each output cell is a weighted mean of its window where the weight of a
//! neighbour is the product of a spatial Gaussian (distance from the centre,
//! `sigma_d`) and a range Gaussian (difference from the centre value,
//! `sigma_r`). Smooths within regions of similar value while leaving sharp
//! transitions intact.

use crate::engine::{map_idx, par_columns, Dims, Edge, Win};
use crate::focal::{check_geom, fill_out, parse_edge};
use extendr_api::prelude::*;

const NAN: f64 = f64::NAN;

#[allow(clippy::too_many_arguments)]
#[inline]
fn bilateral_cell<const NA_AWARE: bool>(
    layer: &[f64],
    d: Dims,
    w: Win,
    spatial: &[f64],
    inv2sr2: f64,
    edge: Edge,
    r: usize,
    c: usize,
    interior: bool,
) -> f64 {
    let nr = d.nr;
    let wr = 2 * w.hr + 1;
    let x0 = layer[c * nr + r];
    if x0.is_nan() {
        return NAN;
    }
    let (mut num, mut den) = (0.0, 0.0);
    if interior {
        for dc in 0..=(2 * w.hc) {
            let base = (c + dc - w.hc) * nr + r - w.hr;
            let scol = &spatial[dc * wr..(dc + 1) * wr];
            for (i, &v) in layer[base..base + wr].iter().enumerate() {
                if v.is_nan() {
                    if NA_AWARE {
                        continue;
                    }
                    return NAN;
                }
                let dv = v - x0;
                let wgt = scol[i] * (-dv * dv * inv2sr2).exp();
                num += wgt * v;
                den += wgt;
            }
        }
    } else {
        for dc in -(w.hc as isize)..=(w.hc as isize) {
            let scol = &spatial[(dc + w.hc as isize) as usize * wr..];
            let cc = map_idx(c as isize + dc, d.nc, edge);
            for dr in -(w.hr as isize)..=(w.hr as isize) {
                let v = match (cc, map_idx(r as isize + dr, nr, edge)) {
                    (Some(cc), Some(rr)) => layer[cc * nr + rr],
                    _ => match edge {
                        Edge::Constant(v) => v,
                        _ => continue,
                    },
                };
                if v.is_nan() {
                    if NA_AWARE {
                        continue;
                    }
                    return NAN;
                }
                let dv = v - x0;
                let wgt =
                    scol[(dr + w.hr as isize) as usize] * (-dv * dv * inv2sr2).exp();
                num += wgt * v;
                den += wgt;
            }
        }
    }
    if den > 0.0 {
        num / den
    } else {
        NAN
    }
}

fn bilateral_impl<const NA_AWARE: bool>(
    x: &[f64],
    d: Dims,
    w: Win,
    sigma_d: f64,
    sigma_r: f64,
    edge: Edge,
    out: &mut [f64],
) {
    let plane = d.nr * d.nc;
    let wr = 2 * w.hr + 1;
    // Spatial weights per window offset, column-major over the window.
    let inv2sd2 = 1.0 / (2.0 * sigma_d * sigma_d);
    let spatial: Vec<f64> = (-(w.hc as isize)..=(w.hc as isize))
        .flat_map(|dc| {
            (-(w.hr as isize)..=(w.hr as isize))
                .map(move |dr| (-((dr * dr + dc * dc) as f64) * inv2sd2).exp())
        })
        .collect();
    debug_assert_eq!(spatial.len(), wr * (2 * w.hc + 1));
    let inv2sr2 = 1.0 / (2.0 * sigma_r * sigma_r);
    par_columns(out, d.nr, |g, out_col| {
        let (l, c) = (g / d.nc, g % d.nc);
        let layer = &x[l * plane..(l + 1) * plane];
        let interior_col = c >= w.hc && c + w.hc < d.nc;
        let r1 = d.nr.saturating_sub(w.hr);
        for (r, o) in out_col.iter_mut().enumerate() {
            let interior = interior_col && r >= w.hr && r < r1;
            *o = bilateral_cell::<NA_AWARE>(layer, d, w, &spatial, inv2sr2, edge, r, c, interior);
        }
    });
}

/// Bilateral filter over a stack of column-major layers.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_bilateral_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    sigma_d: f64,
    sigma_r: f64,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    if !sigma_d.is_finite() || sigma_d <= 0.0 || !sigma_r.is_finite() || sigma_r <= 0.0 {
        throw_r_error("`sigma_d` and `sigma_r` must be positive finite numbers");
    }
    let edge = parse_edge(edge, edge_value);
    fill_out(x.len(), |out| {
        if na_omit {
            bilateral_impl::<true>(x, d, w, sigma_d, sigma_r, edge, out);
        } else {
            bilateral_impl::<false>(x, d, w, sigma_d, sigma_r, edge, out);
        }
    })
}

extendr_module! {
    mod bilateral;
    fn rf_bilateral_rs;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_field_is_invariant() {
        let d = Dims { nr: 5, nc: 5 };
        let x = vec![2.5; 25];
        let mut out = vec![0.0; 25];
        bilateral_impl::<true>(&x, d, Win { hr: 2, hc: 2 }, 1.5, 1.0, Edge::Shrink, &mut out);
        assert!(out.iter().all(|v| (v - 2.5).abs() < 1e-12));
    }

    #[test]
    fn tiny_sigma_r_preserves_a_step_edge() {
        // Two flat regions; with a tiny range sigma each side keeps its value.
        let d = Dims { nr: 4, nc: 6 };
        let x: Vec<f64> = (0..24).map(|i| if i < 12 { 1.0 } else { 10.0 }).collect();
        let mut out = vec![0.0; 24];
        bilateral_impl::<true>(&x, d, Win { hr: 1, hc: 1 }, 1.0, 1e-3, Edge::Shrink, &mut out);
        for (i, v) in out.iter().enumerate() {
            let want = if i < 12 { 1.0 } else { 10.0 };
            assert!((v - want).abs() < 1e-9, "cell {i}: {v}");
        }
    }
}
