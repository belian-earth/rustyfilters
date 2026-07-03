//! SAR speckle filters. All are adaptive local-statistics filters driven by
//! the window mean `m` and population variance `v` around each pixel, with
//! `ci = sd/m` the observed variation coefficient, `cu = 1/sqrt(L)` the
//! speckle variation coefficient for `L`-look intensity data, and
//! `cmax = sqrt(1 + 2/L)` the heterogeneity threshold.
//!
//! Lee, Kuan, enhanced Lee and Gamma-MAP need only the local moments, so they
//! run on the separable O(1)-per-pixel moments path. Frost is two-stage
//! (moments to build the per-pixel damping, then a weighted gather); Lee
//! sigma needs the window values themselves, so both use the scan path.
//!
//! Speckle filters return NA when the centre pixel is NA: the estimators
//! reconstruct the signal *at* the observed pixel, so there is nothing to
//! filter (unlike the smoothing filters, which still summarise neighbours).

use crate::engine::{map_idx, Dims, Edge, Moments, Win};
use crate::focal::{check_geom, fill_out, parse_edge, run_stat_m};
use crate::threading::{get_num_threads, maybe_par};
use extendr_api::prelude::*;
use rayon::prelude::*;

const NAN: f64 = f64::NAN;

/// Guard shared by the single-pass filters: NA centre or empty window -> NA.
/// Returns the window mean and coefficient-of-variation guard state.
#[inline]
fn local_stats(x0: f64, m: &Moments) -> Option<(f64, f64)> {
    if x0.is_nan() || m.n == 0 || m.sum.is_nan() {
        return None;
    }
    let mean = m.mean();
    Some((mean, m.var_pop()))
}

/// Lee (1980): minimum mean-square-error linear estimator.
/// `W = v / (v + (m * cu)^2)`, `out = m + W * (x - m)`.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_lee_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    looks: f64,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    let edge = parse_edge(edge, edge_value);
    let cu2 = 1.0 / looks;
    fill_out(x.len(), |out| {
        run_stat_m(x, d, w, edge, na_omit, |x0, m| {
            match local_stats(x0, &m) {
                None => NAN,
                Some((mean, v)) => {
                    let noise = mean * mean * cu2;
                    if v + noise <= 0.0 {
                        return mean;
                    }
                    let wgt = v / (v + noise);
                    mean + wgt * (x0 - mean)
                }
            }
        }, out);
    })
}

/// Kuan (1985): `W = (1 - cu^2 / ci^2) / (1 + cu^2)` clamped to [0, 1].
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_kuan_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    looks: f64,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    let edge = parse_edge(edge, edge_value);
    let cu2 = 1.0 / looks;
    fill_out(x.len(), |out| {
        run_stat_m(x, d, w, edge, na_omit, |x0, m| {
            match local_stats(x0, &m) {
                None => NAN,
                Some((mean, v)) => {
                    if mean == 0.0 || v == 0.0 {
                        return mean;
                    }
                    let ci2 = v / (mean * mean);
                    let wgt = ((1.0 - cu2 / ci2) / (1.0 + cu2)).clamp(0.0, 1.0);
                    mean + wgt * (x0 - mean)
                }
            }
        }, out);
    })
}

/// Enhanced Lee (Lopes et al. 1990): three regimes split by ci against cu and
/// cmax, with exponential damping between them.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_enhanced_lee_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    looks: f64,
    damping: f64,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    let edge = parse_edge(edge, edge_value);
    let cu = (1.0 / looks).sqrt();
    let cmax = (1.0 + 2.0 / looks).sqrt();
    fill_out(x.len(), |out| {
        run_stat_m(x, d, w, edge, na_omit, |x0, m| {
            match local_stats(x0, &m) {
                None => NAN,
                Some((mean, v)) => {
                    if mean == 0.0 {
                        return mean;
                    }
                    let ci = v.sqrt() / mean;
                    if ci <= cu {
                        mean
                    } else if ci < cmax {
                        let wgt = (-damping * (ci - cu) / (cmax - ci)).exp();
                        mean * wgt + x0 * (1.0 - wgt)
                    } else {
                        x0
                    }
                }
            }
        }, out);
    })
}

/// Gamma-MAP (Lopes et al. 1990): maximum a posteriori under Gamma-
/// distributed signal and speckle.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_gamma_map_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    looks: f64,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    let edge = parse_edge(edge, edge_value);
    let cu2 = 1.0 / looks;
    let cu = cu2.sqrt();
    let cmax = (1.0 + 2.0 / looks).sqrt();
    fill_out(x.len(), |out| {
        run_stat_m(x, d, w, edge, na_omit, |x0, m| {
            match local_stats(x0, &m) {
                None => NAN,
                Some((mean, v)) => {
                    if mean == 0.0 {
                        return mean;
                    }
                    let ci2 = v / (mean * mean);
                    let ci = ci2.sqrt();
                    if ci <= cu {
                        mean
                    } else if ci >= cmax {
                        x0
                    } else {
                        let alpha = (1.0 + cu2) / (ci2 - cu2);
                        let b = alpha - looks - 1.0;
                        let dsc = mean * mean * b * b + 4.0 * alpha * looks * mean * x0;
                        (b * mean + dsc.max(0.0).sqrt()) / (2.0 * alpha)
                    }
                }
            }
        }, out);
    })
}

/// Frost (1982): exponentially damped kernel whose decay adapts to the local
/// variation coefficient. Two stages: the moments path computes the per-pixel
/// damping factor `B = damping * ci^2`, then a scan applies the kernel
/// `exp(-B * dist)` over the window.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_frost_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    damping: f64,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    let edge = parse_edge(edge, edge_value);
    // Stage 1: per-pixel damping factor B = damping * ci^2 (NaN when the
    // centre is NA or, in propagate mode, when the window holds any NA).
    let mut bfac = vec![0.0; x.len()];
    run_stat_m(x, d, w, edge, na_omit, |x0, m| {
        match local_stats(x0, &m) {
            None => NAN,
            Some((mean, v)) => {
                if mean == 0.0 {
                    0.0
                } else {
                    damping * (v / (mean * mean))
                }
            }
        }
    }, &mut bfac);
    // Stage 2: weighted gather with the damped-exponential kernel.
    fill_out(x.len(), |out| {
        if na_omit {
            frost_apply::<true>(x, &bfac, d, w, edge, out);
        } else {
            frost_apply::<false>(x, &bfac, d, w, edge, out);
        }
    })
}

/// Frost stage 2: for each cell, a weighted mean of its window with weights
/// `exp(-B * distance)`. Needs value-offset alignment, so it has its own
/// gather rather than reusing the generic scan.
fn frost_apply<const NA_AWARE: bool>(
    x: &[f64],
    bfac: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    out: &mut [f64],
) {
    let nr = d.nr;
    let plane = nr * d.nc;
    // Distance from the centre per window offset, column-major over the
    // window: dist[(dc + hc) * wr + (dr + hr)].
    let wr = 2 * w.hr + 1;
    let dist: Vec<f64> = (-(w.hc as isize)..=(w.hc as isize))
        .flat_map(|dc| {
            (-(w.hr as isize)..=(w.hr as isize))
                .map(move |dr| ((dr * dr + dc * dc) as f64).sqrt())
        })
        .collect();
    let cell = |layer: &[f64], b: f64, r: usize, c: usize, interior: bool| -> f64 {
        let x0 = layer[c * nr + r];
        if x0.is_nan() || b.is_nan() {
            return f64::NAN;
        }
        let (mut num, mut den) = (0.0, 0.0);
        if interior {
            for dc in 0..=(2 * w.hc) {
                let base = (c + dc - w.hc) * nr + r - w.hr;
                let dcol = &dist[dc * wr..(dc + 1) * wr];
                for (i, &v) in layer[base..base + wr].iter().enumerate() {
                    if v.is_nan() {
                        if NA_AWARE {
                            continue;
                        }
                        return f64::NAN;
                    }
                    let wgt = (-b * dcol[i]).exp();
                    num += wgt * v;
                    den += wgt;
                }
            }
        } else {
            for dc in -(w.hc as isize)..=(w.hc as isize) {
                let dcol = &dist[(dc + w.hc as isize) as usize * wr..];
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
                        return f64::NAN;
                    }
                    let wgt = (-b * dcol[(dr + w.hr as isize) as usize]).exp();
                    num += wgt * v;
                    den += wgt;
                }
            }
        }
        if den > 0.0 {
            num / den
        } else {
            f64::NAN
        }
    };
    let column = |g: usize, out_col: &mut [f64]| {
        let (l, c) = (g / d.nc, g % d.nc);
        let layer = &x[l * plane..(l + 1) * plane];
        let bcol = &bfac[l * plane + c * nr..l * plane + (c + 1) * nr];
        let interior_col = c >= w.hc && c + w.hc < d.nc;
        let r1 = nr.saturating_sub(w.hr);
        for (r, o) in out_col.iter_mut().enumerate() {
            let interior = interior_col && r >= w.hr && r < r1;
            *o = cell(layer, bcol[r], r, c, interior);
        }
    };
    if get_num_threads() <= 1 {
        for (g, out_col) in out.chunks_mut(nr).enumerate() {
            column(g, out_col);
        }
    } else {
        maybe_par(|| {
            out.par_chunks_mut(nr)
                .enumerate()
                .for_each(|(g, out_col)| column(g, out_col));
        });
    }
}

/// Lee sigma (1983): mean of the window pixels inside the two-sigma bounds
/// `x * (1 +/- k * cu)`; if fewer than `min_count` pixels qualify, fall back
/// to the full window mean.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_lee_sigma_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    looks: f64,
    k: f64,
    min_count: i32,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    let edge = parse_edge(edge, edge_value);
    let cu = (1.0 / looks).sqrt();
    let lo_f = 1.0 - k * cu;
    let hi_f = 1.0 + k * cu;
    let min_count = min_count.max(1) as usize;
    fill_out(x.len(), |out| {
        crate::focal::run_stat(x, d, w, edge, na_omit, |x0, v, m| {
            if x0.is_nan() || m.n == 0 || m.sum.is_nan() {
                return NAN;
            }
            // Intensity data is non-negative; order the bounds anyway so a
            // negative centre still gets a well-formed interval.
            let (lo, hi) = if x0 >= 0.0 {
                (x0 * lo_f, x0 * hi_f)
            } else {
                (x0 * hi_f, x0 * lo_f)
            };
            let (mut s, mut n) = (0.0, 0usize);
            for &val in v.iter() {
                if val >= lo && val <= hi {
                    s += val;
                    n += 1;
                }
            }
            if n < min_count {
                m.mean()
            } else {
                s / n as f64
            }
        }, out);
    })
}

extendr_module! {
    mod speckle;
    fn rf_lee_rs;
    fn rf_kuan_rs;
    fn rf_enhanced_lee_rs;
    fn rf_gamma_map_rs;
    fn rf_frost_rs;
    fn rf_lee_sigma_rs;
}
