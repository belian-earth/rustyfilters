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

use crate::engine::{map_idx, par_columns, Dims, Edge, Moments, Win};
use crate::focal::{check_geom, fill_out, parse_edge, run_stat_m};
use extendr_api::prelude::*;

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
        run_stat_m(
            x,
            d,
            w,
            edge,
            na_omit,
            |x0, m| match local_stats(x0, &m) {
                None => NAN,
                Some((mean, v)) => {
                    let noise = mean * mean * cu2;
                    if v + noise <= 0.0 {
                        return mean;
                    }
                    let wgt = v / (v + noise);
                    mean + wgt * (x0 - mean)
                }
            },
            out,
        );
    })
}

/// Kuan (1985): `W = (1 - cu^2 / ci^2) / (1 + cu^2)` clamped to `[0, 1]`.
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
        run_stat_m(
            x,
            d,
            w,
            edge,
            na_omit,
            |x0, m| match local_stats(x0, &m) {
                None => NAN,
                Some((mean, v)) => {
                    if mean == 0.0 || v == 0.0 {
                        return mean;
                    }
                    let ci2 = v / (mean * mean);
                    let wgt = ((1.0 - cu2 / ci2) / (1.0 + cu2)).clamp(0.0, 1.0);
                    mean + wgt * (x0 - mean)
                }
            },
            out,
        );
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
        run_stat_m(
            x,
            d,
            w,
            edge,
            na_omit,
            |x0, m| match local_stats(x0, &m) {
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
            },
            out,
        );
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
        run_stat_m(
            x,
            d,
            w,
            edge,
            na_omit,
            |x0, m| match local_stats(x0, &m) {
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
            },
            out,
        );
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
    run_stat_m(
        x,
        d,
        w,
        edge,
        na_omit,
        |x0, m| match local_stats(x0, &m) {
            None => NAN,
            Some((mean, v)) => {
                if mean == 0.0 {
                    0.0
                } else {
                    damping * (v / (mean * mean))
                }
            }
        },
        &mut bfac,
    );
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
            (-(w.hr as isize)..=(w.hr as isize)).map(move |dr| ((dr * dr + dc * dc) as f64).sqrt())
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
    par_columns(out, nr, |g, out_col| {
        let (l, c) = (g / d.nc, g % d.nc);
        let layer = &x[l * plane..(l + 1) * plane];
        let bcol = &bfac[l * plane + c * nr..l * plane + (c + 1) * nr];
        let interior_col = c >= w.hc && c + w.hc < d.nc;
        let r1 = nr.saturating_sub(w.hr);
        for (r, o) in out_col.iter_mut().enumerate() {
            let interior = interior_col && r >= w.hr && r < r1;
            *o = cell(layer, bcol[r], r, c, interior);
        }
    });
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
        crate::focal::run_stat(
            x,
            d,
            w,
            edge,
            na_omit,
            |x0, v, m| {
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
            },
            out,
        );
    })
}

// ---------------------------------------------------------------------------
// Improved Lee sigma (Lee et al. 2009).
//
// Sigma range bounds (i1, i2) and the revised speckle variation coefficient
// eta_v' for intensity data, from Table I/II of Lee, Wen, Ainsworth, Chen &
// Chen (2009), "Improved sigma filter for speckle filtering of SAR imagery",
// IEEE TGRS 47(1). Rows: looks 1-4; columns: sigma 0.5-0.9.
const SIGMA_2009: [[(f64, f64, f64); 5]; 4] = [
    [
        (0.436, 1.920, 0.4057),
        (0.343, 2.210, 0.4954),
        (0.254, 2.582, 0.5911),
        (0.168, 3.094, 0.6966),
        (0.084, 3.941, 0.8191),
    ],
    [
        (0.582, 1.584, 0.2763),
        (0.501, 1.755, 0.3388),
        (0.418, 1.972, 0.4062),
        (0.327, 2.260, 0.4810),
        (0.221, 2.744, 0.5699),
    ],
    [
        (0.652, 1.458, 0.2222),
        (0.580, 1.586, 0.2736),
        (0.505, 1.751, 0.3280),
        (0.419, 1.965, 0.3892),
        (0.313, 2.320, 0.4624),
    ],
    [
        (0.694, 1.385, 0.1921),
        (0.630, 1.495, 0.2348),
        (0.560, 1.627, 0.2825),
        (0.480, 1.804, 0.3354),
        (0.378, 2.094, 0.3991),
    ],
];

/// Bright pixels within a target window needed to call a cluster a point
/// target (Lee et al. 2009 / SNAP convention).
const TARGET_CLUSTER_SIZE: usize = 5;

/// MMSE estimate: `(1 - b) * mean + b * x0` with `b = varX / varY`,
/// `varX = (varY - mean^2 * eta2) / (1 + eta2)` clamped at 0. Sample
/// variance, matching the 2009 filter's reference implementation.
fn mmse(vals: &[f64], x0: f64, eta2: f64) -> f64 {
    let n = vals.len();
    let mean = vals.iter().sum::<f64>() / n as f64;
    if n < 2 {
        return mean;
    }
    let var_y = vals.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1) as f64;
    if var_y == 0.0 {
        return mean;
    }
    let var_x = ((var_y - mean * mean * eta2) / (1.0 + eta2)).max(0.0);
    let b = var_x / var_y;
    (1.0 - b) * mean + b * x0
}

/// Gather the valid values of a window around (r, c), optionally keeping only
/// those inside `range`. Returns None in propagate mode when the window holds
/// a NaN.
#[allow(clippy::too_many_arguments)]
fn gather_valid<const NA_AWARE: bool>(
    layer: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    r: usize,
    c: usize,
    range: Option<(f64, f64)>,
    scratch: &mut Vec<f64>,
) -> Option<()> {
    scratch.clear();
    for dc in -(w.hc as isize)..=(w.hc as isize) {
        let cc = map_idx(c as isize + dc, d.nc, edge);
        for dr in -(w.hr as isize)..=(w.hr as isize) {
            let v = match (cc, map_idx(r as isize + dr, d.nr, edge)) {
                (Some(cc), Some(rr)) => layer[cc * d.nr + rr],
                _ => match edge {
                    Edge::Constant(v) => v,
                    _ => continue,
                },
            };
            if v.is_nan() {
                if NA_AWARE {
                    continue;
                }
                return None;
            }
            match range {
                Some((lo, hi)) if v < lo || v > hi => {}
                _ => scratch.push(v),
            }
        }
    }
    Some(())
}

/// 98th percentile of the valid cells of one layer.
fn z98_of(layer: &[f64]) -> f64 {
    let mut vals: Vec<f64> = layer.iter().copied().filter(|v| !v.is_nan()).collect();
    if vals.is_empty() {
        return f64::INFINITY;
    }
    let idx = (((vals.len() as f64) * 0.98) as usize)
        .saturating_sub(1)
        .min(vals.len() - 1);
    let (_, &mut v, _) = vals.select_nth_unstable_by(idx, |a, b| a.total_cmp(b));
    v
}

#[allow(clippy::too_many_arguments)]
fn lee_sigma_improved_layer<const NA_AWARE: bool>(
    layer: &[f64],
    d: Dims,
    w: Win,
    tw: Win,
    i1: f64,
    i2: f64,
    eta_v2: f64,
    eta_vp2: f64,
    edge: Edge,
    out: &mut [f64],
) {
    let nr = d.nr;
    let z98 = z98_of(layer);
    // Pass 1: detectors -- bright pixels whose target window holds a bright
    // cluster. Pass 2: point targets -- bright pixels near a detector. This
    // is an order-independent variant of the 2009 paper's sequential
    // marking; detected clusters are preserved untouched either way.
    let mut detector = vec![0.0; layer.len()];
    par_columns(&mut detector, nr, |c, det_col| {
        for (r, o) in det_col.iter_mut().enumerate() {
            let v = layer[c * nr + r];
            if v.is_nan() || v <= z98 {
                continue;
            }
            let mut cluster = 0usize;
            for dc in -(tw.hc as isize)..=(tw.hc as isize) {
                let cc = map_idx(c as isize + dc, d.nc, edge);
                for dr in -(tw.hr as isize)..=(tw.hr as isize) {
                    if let (Some(cc), Some(rr)) = (cc, map_idx(r as isize + dr, nr, edge)) {
                        let q = layer[cc * nr + rr];
                        if !q.is_nan() && q > z98 {
                            cluster += 1;
                        }
                    }
                }
            }
            if cluster > TARGET_CLUSTER_SIZE {
                *o = 1.0;
            }
        }
    });
    let mut target = vec![0.0; layer.len()];
    par_columns(&mut target, nr, |c, tgt_col| {
        for (r, o) in tgt_col.iter_mut().enumerate() {
            let v = layer[c * nr + r];
            if v.is_nan() || v <= z98 {
                continue;
            }
            'search: for dc in -(tw.hc as isize)..=(tw.hc as isize) {
                let cc = map_idx(c as isize + dc, d.nc, edge);
                for dr in -(tw.hr as isize)..=(tw.hr as isize) {
                    if let (Some(cc), Some(rr)) = (cc, map_idx(r as isize + dr, nr, edge)) {
                        if detector[cc * nr + rr] == 1.0 {
                            *o = 1.0;
                            break 'search;
                        }
                    }
                }
            }
        }
    });
    // Pass 3: filter.
    par_columns(out, nr, |c, out_col| {
        let mut scratch = Vec::with_capacity((2 * w.hr + 1) * (2 * w.hc + 1));
        for (r, o) in out_col.iter_mut().enumerate() {
            let x0 = layer[c * nr + r];
            if x0.is_nan() {
                *o = f64::NAN;
                continue;
            }
            if target[c * nr + r] == 1.0 {
                *o = x0;
                continue;
            }
            // A priori mean from the target window, then MMSE over the
            // filter-window pixels inside the sigma range.
            let est = match gather_valid::<NA_AWARE>(layer, d, tw, edge, r, c, None, &mut scratch)
            {
                None => {
                    *o = f64::NAN;
                    continue;
                }
                Some(()) if scratch.is_empty() => x0,
                Some(()) => mmse(&scratch, x0, eta_v2),
            };
            match gather_valid::<NA_AWARE>(
                layer,
                d,
                w,
                edge,
                r,
                c,
                Some((est * i1, est * i2)),
                &mut scratch,
            ) {
                None => *o = f64::NAN,
                Some(()) if scratch.is_empty() => *o = x0,
                Some(()) => *o = mmse(&scratch, x0, eta_vp2),
            }
        }
    });
}

/// Improved Lee sigma (Lee et al. 2009) over a stack of column-major layers.
/// `sigma_idx` indexes 0.5-0.9; `looks` must be 1-4; `twr`/`twc` give the
/// target window.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_lee_sigma_improved_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    looks: i32,
    sigma_idx: i32,
    twr: i32,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    if !(1..=4).contains(&looks) {
        throw_r_error("`looks` must be 1, 2, 3 or 4 for the improved Lee sigma filter");
    }
    if !(0..=4).contains(&sigma_idx) {
        throw_r_error("internal error: bad sigma index");
    }
    if twr < 1 || twr % 2 == 0 {
        throw_r_error("`target_window` must be an odd positive integer");
    }
    let tw = Win {
        hr: (twr as usize - 1) / 2,
        hc: (twr as usize - 1) / 2,
    };
    let (i1, i2, eta_vp) = SIGMA_2009[(looks - 1) as usize][sigma_idx as usize];
    let eta_v2 = 1.0 / looks as f64;
    let eta_vp2 = eta_vp * eta_vp;
    let edge = parse_edge(edge, edge_value);
    let plane = d.nr * d.nc;
    fill_out(x.len(), |out| {
        for (layer, out_layer) in x.chunks(plane).zip(out.chunks_mut(plane)) {
            if na_omit {
                lee_sigma_improved_layer::<true>(
                    layer, d, w, tw, i1, i2, eta_v2, eta_vp2, edge, out_layer,
                );
            } else {
                lee_sigma_improved_layer::<false>(
                    layer, d, w, tw, i1, i2, eta_v2, eta_vp2, edge, out_layer,
                );
            }
        }
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
    fn rf_lee_sigma_improved_rs;
}
