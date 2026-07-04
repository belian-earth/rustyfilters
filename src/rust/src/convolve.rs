//! Generic kernel filtering: the `rf_convolve()` backend (and, via fixed
//! kernels on the R side, `rf_sobel()` and `rf_laplacian()`).
//!
//! The kernel is applied as-is (cross-correlation): for a `kr x kc` kernel K
//! with half sizes (hr, hc), `out(r, c) = sum K[dr + hr, dc + hc] *
//! x(r + dr, c + dc)`. Flip the kernel beforehand for a true convolution;
//! the distinction only matters for asymmetric kernels.

use crate::engine::{map_idx, par_columns, Dims, Edge, Win};
use crate::focal::{check_geom, fill_out, parse_edge};
use extendr_api::prelude::*;

const NAN: f64 = f64::NAN;

/// Weighted gather for one cell. Under `NA_AWARE` missing cells contribute
/// nothing; otherwise any missing cell makes the result NA. `normalize`
/// divides by the sum of the weights actually used (normalized convolution).
#[allow(clippy::too_many_arguments)]
#[inline]
fn conv_cell<const NA_AWARE: bool>(
    layer: &[f64],
    d: Dims,
    w: Win,
    kernel: &[f64],
    edge: Edge,
    normalize: bool,
    r: usize,
    c: usize,
    interior: bool,
) -> f64 {
    let nr = d.nr;
    let kr = 2 * w.hr + 1;
    let (mut num, mut den) = (0.0, 0.0);
    let mut n = 0usize;
    if interior {
        for dc in 0..=(2 * w.hc) {
            let base = (c + dc - w.hc) * nr + r - w.hr;
            let kcol = &kernel[dc * kr..(dc + 1) * kr];
            for (i, &v) in layer[base..base + kr].iter().enumerate() {
                if v.is_nan() {
                    if NA_AWARE {
                        continue;
                    }
                    return NAN;
                }
                num += kcol[i] * v;
                den += kcol[i];
                n += 1;
            }
        }
    } else {
        for dc in -(w.hc as isize)..=(w.hc as isize) {
            let kcol = &kernel[(dc + w.hc as isize) as usize * kr..];
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
                let kw = kcol[(dr + w.hr as isize) as usize];
                num += kw * v;
                den += kw;
                n += 1;
            }
        }
    }
    if n == 0 {
        return NAN;
    }
    if normalize {
        if den.abs() > 1e-12 {
            num / den
        } else {
            NAN
        }
    } else {
        num
    }
}

fn convolve_impl<const NA_AWARE: bool>(
    x: &[f64],
    d: Dims,
    w: Win,
    kernel: &[f64],
    edge: Edge,
    normalize: bool,
    out: &mut [f64],
) {
    let plane = d.nr * d.nc;
    par_columns(out, d.nr, |g, out_col| {
        let (l, c) = (g / d.nc, g % d.nc);
        let layer = &x[l * plane..(l + 1) * plane];
        let interior_col = c >= w.hc && c + w.hc < d.nc;
        let r1 = d.nr.saturating_sub(w.hr);
        for (r, o) in out_col.iter_mut().enumerate() {
            let interior = interior_col && r >= w.hr && r < r1;
            *o = conv_cell::<NA_AWARE>(layer, d, w, kernel, edge, normalize, r, c, interior);
        }
    });
}

/// Apply a kernel (column-major, `wr x wc`) over a stack of layers.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_convolve_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    kernel: &[f64],
    normalize: bool,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    if kernel.len() != (wr as usize) * (wc as usize) {
        throw_r_error("kernel length does not match its dimensions");
    }
    let edge = parse_edge(edge, edge_value);
    fill_out(x.len(), |out| {
        if na_omit {
            convolve_impl::<true>(x, d, w, kernel, edge, normalize, out);
        } else {
            convolve_impl::<false>(x, d, w, kernel, edge, normalize, out);
        }
    })
}

extendr_module! {
    mod convolve;
    fn rf_convolve_rs;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_kernel_returns_input() {
        let d = Dims { nr: 4, nc: 5 };
        let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let k = [0., 0., 0., 0., 1., 0., 0., 0., 0.];
        let mut out = vec![0.0; 20];
        convolve_impl::<true>(&x, d, Win { hr: 1, hc: 1 }, &k, Edge::Shrink, false, &mut out);
        assert_eq!(out, x);
    }

    #[test]
    fn box_kernel_with_normalize_is_the_mean() {
        let d = Dims { nr: 3, nc: 3 };
        let x = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let k = [1.0; 9];
        let mut out = vec![0.0; 9];
        convolve_impl::<true>(&x, d, Win { hr: 1, hc: 1 }, &k, Edge::Shrink, true, &mut out);
        assert!((out[4] - 5.0).abs() < 1e-12);
        assert!((out[0] - 3.0).abs() < 1e-12); // corner: mean of 4 cells
    }

    #[test]
    fn fast_path_poisons_na_windows() {
        let d = Dims { nr: 3, nc: 3 };
        let x = [1., f64::NAN, 3., 4., 5., 6., 7., 8., 9.];
        let k = [1.0; 9];
        let mut out = vec![0.0; 9];
        convolve_impl::<false>(&x, d, Win { hr: 1, hc: 1 }, &k, Edge::Shrink, false, &mut out);
        assert!(out[4].is_nan());
        assert!(!out[8].is_nan());
    }
}
