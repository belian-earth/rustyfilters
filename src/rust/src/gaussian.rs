//! Separable Gaussian blur: the `rf_gaussian()` backend.
//!
//! NA-aware mode is normalized convolution: convolve `x` with NaNs zeroed
//! and the validity mask through the same separable kernel, then divide.
//! Numerator and denominator are each separable, so the two-pass result
//! equals the full 2-D normalized convolution even with arbitrary NA holes.
//! The fast path convolves the data alone; `shrink` edges are renormalized
//! by the analytic marginal kernel sums (the mask of pure border truncation
//! is separable, so no mask pass is needed).

use crate::engine::{map_idx, Dims, Edge};
use crate::focal::{check_geom, fill_out, parse_edge};
use crate::threading::{get_num_threads, maybe_par};
use extendr_api::prelude::*;
use rayon::prelude::*;

const NAN: f64 = f64::NAN;

/// Normalized Gaussian weights for an odd window length.
fn gauss_kernel(w: usize, sigma: f64) -> Vec<f64> {
    let h = (w / 2) as isize;
    let mut k: Vec<f64> = (-h..=h)
        .map(|i| (-((i * i) as f64) / (2.0 * sigma * sigma)).exp())
        .collect();
    let s: f64 = k.iter().sum();
    k.iter_mut().for_each(|v| *v /= s);
    k
}

/// One cell's contribution to the (numerator, weight) pair.
#[inline]
fn acc_weighted<const NA_AWARE: bool>(v: f64, kw: f64, num: &mut f64, den: &mut f64) {
    if NA_AWARE && v.is_nan() {
        return;
    }
    *num += kw * v;
    *den += kw;
}

/// Vertical pass over one column: `tx[r] = sum k[i] * x[r+i-h]`, `tm[r]` the
/// matching kernel-weight sum over valid cells (NA-aware mode only).
fn conv_col<const NA_AWARE: bool>(
    col: &[f64],
    nr: usize,
    k: &[f64],
    edge: Edge,
    tx: &mut [f64],
    tm: &mut [f64],
) {
    let h = k.len() / 2;
    let direct = |r: usize, tx: &mut [f64], tm: &mut [f64]| {
        let (mut num, mut den) = (0.0, 0.0);
        for (i, &kw) in k.iter().enumerate() {
            match map_idx(r as isize + i as isize - h as isize, nr, edge) {
                Some(rr) => acc_weighted::<NA_AWARE>(col[rr], kw, &mut num, &mut den),
                None => {
                    if let Edge::Constant(v) = edge {
                        acc_weighted::<NA_AWARE>(v, kw, &mut num, &mut den);
                    }
                }
            }
        }
        tx[r] = num;
        tm[r] = den;
    };
    let r1 = nr.saturating_sub(h);
    if h >= r1 {
        for r in 0..nr {
            direct(r, tx, tm);
        }
        return;
    }
    for r in 0..h {
        direct(r, tx, tm);
    }
    for r in h..r1 {
        let (mut num, mut den) = (0.0, 0.0);
        for (i, &kw) in k.iter().enumerate() {
            acc_weighted::<NA_AWARE>(col[r - h + i], kw, &mut num, &mut den);
        }
        tx[r] = num;
        tm[r] = den;
    }
    for r in r1..nr {
        direct(r, tx, tm);
    }
}

/// Horizontal pass for one output column plus the final division.
#[allow(clippy::too_many_arguments)]
fn finish_col<const NA_AWARE: bool>(
    tx: &[f64],
    tm: &[f64],
    d: Dims,
    k: &[f64],
    edge: Edge,
    c: usize,
    vmarg: &[f64],
    hmarg: &[f64],
    out_col: &mut [f64],
) {
    let h = k.len() / 2;
    let nr = d.nr;
    let mut num = vec![0.0; nr];
    let mut den = vec![0.0; nr];
    for (j, &kw) in k.iter().enumerate() {
        match map_idx(c as isize + j as isize - h as isize, d.nc, edge) {
            Some(cc) => {
                let txc = &tx[cc * nr..(cc + 1) * nr];
                if NA_AWARE {
                    let tmc = &tm[cc * nr..(cc + 1) * nr];
                    for r in 0..nr {
                        num[r] += kw * txc[r];
                        den[r] += kw * tmc[r];
                    }
                } else {
                    for r in 0..nr {
                        num[r] += kw * txc[r];
                    }
                }
            }
            None => {
                // Whole kernel column outside the grid.
                if let Edge::Constant(v) = edge {
                    if NA_AWARE && v.is_nan() {
                        continue;
                    }
                    // A constant column convolves vertically to (v, 1).
                    num.iter_mut().for_each(|n| *n += kw * v);
                    if NA_AWARE {
                        den.iter_mut().for_each(|d| *d += kw);
                    }
                }
            }
        }
    }
    if NA_AWARE {
        for r in 0..nr {
            out_col[r] = if den[r] > 1e-12 { num[r] / den[r] } else { NAN };
        }
    } else if matches!(edge, Edge::Shrink) {
        // Border truncation renormalizes by the separable marginal sums.
        let hm = hmarg[c];
        for r in 0..nr {
            out_col[r] = num[r] / (vmarg[r] * hm);
        }
    } else {
        out_col.copy_from_slice(&num);
    }
}

/// Kernel-weight sum actually inside the axis, per position (all 1 away from
/// the borders since kernels are normalized).
fn marginal_sums(n: usize, k: &[f64]) -> Vec<f64> {
    let h = k.len() / 2;
    (0..n)
        .map(|p| {
            k.iter()
                .enumerate()
                .filter(|(i, _)| {
                    let q = p as isize + *i as isize - h as isize;
                    q >= 0 && (q as usize) < n
                })
                .map(|(_, &kw)| kw)
                .sum()
        })
        .collect()
}

fn gaussian_impl<const NA_AWARE: bool>(
    x: &[f64],
    d: Dims,
    kr: &[f64],
    kc: &[f64],
    edge: Edge,
    out: &mut [f64],
) {
    let plane = d.nr * d.nc;
    let seq = get_num_threads() <= 1;
    let vmarg = marginal_sums(d.nr, kr);
    let hmarg = marginal_sums(d.nc, kc);
    let mut tx = vec![0.0; plane];
    let mut tm = vec![0.0; if NA_AWARE { plane } else { 0 }];
    // In the fast path the mask pass is skipped; give conv_col a scratch row
    // so the write target exists without a branch in the hot loop.
    let mut tm_dummy = vec![0.0; if NA_AWARE { 0 } else { d.nr }];
    for (layer_x, layer_out) in x.chunks(plane).zip(out.chunks_mut(plane)) {
        // Pass 1: vertical, one task per column.
        if NA_AWARE {
            let cols = tx.chunks_mut(d.nr).zip(tm.chunks_mut(d.nr)).enumerate();
            if seq {
                for (c, (txc, tmc)) in cols {
                    conv_col::<NA_AWARE>(
                        &layer_x[c * d.nr..(c + 1) * d.nr],
                        d.nr,
                        kr,
                        edge,
                        txc,
                        tmc,
                    );
                }
            } else {
                maybe_par(|| {
                    tx.par_chunks_mut(d.nr)
                        .zip(tm.par_chunks_mut(d.nr))
                        .enumerate()
                        .for_each(|(c, (txc, tmc))| {
                            conv_col::<NA_AWARE>(
                                &layer_x[c * d.nr..(c + 1) * d.nr],
                                d.nr,
                                kr,
                                edge,
                                txc,
                                tmc,
                            );
                        });
                });
            }
        } else if seq {
            for (c, txc) in tx.chunks_mut(d.nr).enumerate() {
                conv_col::<NA_AWARE>(
                    &layer_x[c * d.nr..(c + 1) * d.nr],
                    d.nr,
                    kr,
                    edge,
                    txc,
                    &mut tm_dummy,
                );
            }
        } else {
            maybe_par(|| {
                tx.par_chunks_mut(d.nr).enumerate().for_each_init(
                    || vec![0.0; d.nr],
                    |tmd, (c, txc)| {
                        conv_col::<NA_AWARE>(
                            &layer_x[c * d.nr..(c + 1) * d.nr],
                            d.nr,
                            kr,
                            edge,
                            txc,
                            tmd,
                        );
                    },
                );
            });
        }
        // Pass 2: horizontal + divide, one task per output column.
        let finish = |c: usize, out_col: &mut [f64]| {
            finish_col::<NA_AWARE>(&tx, &tm, d, kc, edge, c, &vmarg, &hmarg, out_col);
        };
        if seq {
            for (c, out_col) in layer_out.chunks_mut(d.nr).enumerate() {
                finish(c, out_col);
            }
        } else {
            maybe_par(|| {
                layer_out
                    .par_chunks_mut(d.nr)
                    .enumerate()
                    .for_each(|(c, out_col)| finish(c, out_col));
            });
        }
    }
}

/// Separable Gaussian blur over a stack of column-major layers.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_gaussian_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    sigma: f64,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    if !sigma.is_finite() || sigma <= 0.0 {
        throw_r_error("`sigma` must be a positive finite number");
    }
    let edge = parse_edge(edge, edge_value);
    let kr = gauss_kernel(2 * w.hr + 1, sigma);
    let kc = gauss_kernel(2 * w.hc + 1, sigma);
    fill_out(x.len(), |out_f| {
        if na_omit {
            gaussian_impl::<true>(x, d, &kr, &kc, edge, out_f);
        } else {
            gaussian_impl::<false>(x, d, &kr, &kc, edge, out_f);
        }
    })
}

extendr_module! {
    mod gaussian;
    fn rf_gaussian_rs;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kernel_is_normalized_and_symmetric() {
        let k = gauss_kernel(7, 1.5);
        assert!((k.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!((k[0] - k[6]).abs() < 1e-15);
        assert!(k[3] > k[2]);
    }

    #[test]
    fn constant_field_is_invariant() {
        let d = Dims { nr: 6, nc: 5 };
        let x = vec![3.5; 30];
        let k = gauss_kernel(5, 1.0);
        let mut out = vec![0.0; 30];
        for edge in [
            Edge::Shrink,
            Edge::Reflect,
            Edge::Nearest,
            Edge::Constant(3.5),
        ] {
            gaussian_impl::<true>(&x, d, &k, &k, edge, &mut out);
            assert!(out.iter().all(|v| (v - 3.5).abs() < 1e-12), "{edge:?}");
            gaussian_impl::<false>(&x, d, &k, &k, edge, &mut out);
            assert!(out.iter().all(|v| (v - 3.5).abs() < 1e-12), "{edge:?}");
        }
    }

    #[test]
    fn na_hole_is_interpolated_in_omit_mode() {
        let d = Dims { nr: 5, nc: 5 };
        let mut x = vec![2.0; 25];
        x[12] = f64::NAN;
        let k = gauss_kernel(3, 1.0);
        let mut out = vec![0.0; 25];
        gaussian_impl::<true>(&x, d, &k, &k, Edge::Shrink, &mut out);
        assert!((out[12] - 2.0).abs() < 1e-12);
        gaussian_impl::<false>(&x, d, &k, &k, Edge::Shrink, &mut out);
        assert!(out[12].is_nan() && out[6].is_nan() && !out[0].is_nan());
    }
}
