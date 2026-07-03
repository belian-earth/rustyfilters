//! Focal statistics over a moving window: the `rf_focal()` backend.

use crate::engine::{run_scan, Dims, Edge, Moments, Win};
use extendr_api::prelude::*;

/// Parse the edge policy string validated on the R side.
pub(crate) fn parse_edge(edge: &str, edge_value: f64) -> Edge {
    match edge {
        "shrink" => Edge::Shrink,
        "reflect" => Edge::Reflect,
        "nearest" => Edge::Nearest,
        "constant" => Edge::Constant(edge_value),
        _ => throw_r_error(format!("unknown edge policy: {edge}")),
    }
}

/// Validate geometry shared by every filter entry point.
pub(crate) fn check_geom(x_len: usize, nr: i32, nc: i32, nl: i32, wr: i32, wc: i32) -> (Dims, Win) {
    if nr < 1 || nc < 1 || nl < 1 {
        throw_r_error("dimensions must be positive");
    }
    if wr < 1 || wc < 1 || wr % 2 == 0 || wc % 2 == 0 {
        throw_r_error("window sizes must be odd positive integers");
    }
    let (nr, nc, nl) = (nr as usize, nc as usize, nl as usize);
    if nr.checked_mul(nc).and_then(|p| p.checked_mul(nl)) != Some(x_len) {
        throw_r_error("length of `x` does not match its dimensions");
    }
    (
        Dims { nr, nc },
        Win {
            hr: (wr as usize - 1) / 2,
            hc: (wc as usize - 1) / 2,
        },
    )
}

/// Allocate the output R vector, expose it to `fill` as a plain f64 slice,
/// and return it.
///
/// SAFETY: `Rfloat` is `#[repr(transparent)]` over `f64`, so the slice cast
/// is layout-sound. The engine only writes plain memory from worker threads
/// (no R API calls), which is safe on an already-allocated vector.
pub(crate) fn fill_out<G: FnOnce(&mut [f64])>(len: usize, fill: G) -> Doubles {
    let mut out = Doubles::new(len);
    let slice: &mut [Rfloat] = &mut out;
    let out_f = unsafe { std::slice::from_raw_parts_mut(slice.as_mut_ptr() as *mut f64, len) };
    fill(out_f);
    out
}

/// Monomorphise the NA policy branch for a reducer closure.
pub(crate) fn run_stat<F>(
    x: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    na_omit: bool,
    f: F,
    out: &mut [f64],
) where
    F: Fn(f64, &mut Vec<f64>, Moments) -> f64 + Sync,
{
    if na_omit {
        run_scan::<F, true>(x, d, w, edge, &f, out);
    } else {
        run_scan::<F, false>(x, d, w, edge, &f, out);
    }
}

/// Median of a scratch window; even counts average the two middle values.
fn median_of(v: &mut [f64]) -> f64 {
    let n = v.len();
    let mid = n / 2;
    let (_, &mut hi, _) = v.select_nth_unstable_by(mid, |a, b| a.total_cmp(b));
    if n % 2 == 1 {
        hi
    } else {
        let lo = v[..mid]
            .iter()
            .fold(f64::NEG_INFINITY, |a, &b| if b > a { b } else { a });
        (lo + hi) / 2.0
    }
}

/// Most frequent value; ties resolve to the lowest value.
fn mode_of(v: &mut [f64]) -> f64 {
    v.sort_unstable_by(|a, b| a.total_cmp(b));
    let (mut best, mut best_n) = (v[0], 0usize);
    let (mut cur, mut cur_n) = (v[0], 0usize);
    for &x in v.iter() {
        if x == cur {
            cur_n += 1;
        } else {
            if cur_n > best_n {
                (best, best_n) = (cur, cur_n);
            }
            (cur, cur_n) = (x, 1);
        }
    }
    if cur_n > best_n {
        best = cur;
    }
    best
}

const NAN: f64 = f64::NAN;

/// Shared guard for order statistics: empty window (all NA in omit mode) or
/// a poisoned sum (any NA in the fast path) both yield NA.
#[inline]
fn bad_window(v: &[f64], m: &Moments) -> bool {
    v.is_empty() || m.sum.is_nan()
}

pub(crate) fn focal_stat(
    stat: &str,
    x: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    na_omit: bool,
    out: &mut [f64],
) {
    match stat {
        "mean" => run_stat(x, d, w, edge, na_omit, |_, _, m| {
            if m.n == 0 {
                NAN
            } else {
                m.mean()
            }
        }, out),
        "sum" => run_stat(x, d, w, edge, na_omit, |_, _, m| {
            if m.n == 0 {
                NAN
            } else {
                m.sum
            }
        }, out),
        "sd" => run_stat(x, d, w, edge, na_omit, |_, _, m| {
            if m.n < 2 {
                NAN
            } else {
                m.var_samp().sqrt()
            }
        }, out),
        "min" => run_stat(x, d, w, edge, na_omit, |_, v, m| {
            if bad_window(v, &m) {
                NAN
            } else {
                v.iter().fold(f64::INFINITY, |a, &b| if b < a { b } else { a })
            }
        }, out),
        "max" => run_stat(x, d, w, edge, na_omit, |_, v, m| {
            if bad_window(v, &m) {
                NAN
            } else {
                v.iter().fold(f64::NEG_INFINITY, |a, &b| if b > a { b } else { a })
            }
        }, out),
        "range" => run_stat(x, d, w, edge, na_omit, |_, v, m| {
            if bad_window(v, &m) {
                NAN
            } else {
                let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                for &b in v.iter() {
                    if b < lo {
                        lo = b;
                    }
                    if b > hi {
                        hi = b;
                    }
                }
                hi - lo
            }
        }, out),
        "median" => run_stat(x, d, w, edge, na_omit, |_, v, m| {
            if bad_window(v, &m) {
                NAN
            } else {
                median_of(v)
            }
        }, out),
        "mode" => run_stat(x, d, w, edge, na_omit, |_, v, m| {
            if bad_window(v, &m) {
                NAN
            } else {
                mode_of(v)
            }
        }, out),
        _ => throw_r_error(format!("unknown stat: {stat}")),
    }
}

/// Focal statistic over a stack of column-major layers.
/// @noRd
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rf_focal_rs(
    x: &[f64],
    nr: i32,
    nc: i32,
    nl: i32,
    wr: i32,
    wc: i32,
    stat: &str,
    edge: &str,
    edge_value: f64,
    na_omit: bool,
) -> Doubles {
    let (d, w) = check_geom(x.len(), nr, nc, nl, wr, wc);
    let edge = parse_edge(edge, edge_value);
    fill_out(x.len(), |out_f| {
        focal_stat(stat, x, d, w, edge, na_omit, out_f)
    })
}

extendr_module! {
    mod focal;
    fn rf_focal_rs;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_even_and_odd() {
        let mut v = [3., 1., 2.];
        assert_eq!(median_of(&mut v), 2.0);
        let mut v = [4., 1., 3., 2.];
        assert_eq!(median_of(&mut v), 2.5);
    }

    #[test]
    fn mode_ties_take_lowest() {
        let mut v = [2., 1., 2., 1., 3.];
        assert_eq!(mode_of(&mut v), 1.0);
        let mut v = [5., 5., 1.];
        assert_eq!(mode_of(&mut v), 5.0);
    }
}
