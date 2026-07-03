//! Core moving-window machinery shared by every filter. Pure Rust (no
//! extendr types) so it is unit-testable with plain `cargo test`.
//!
//! Data is column-major (R layout): a layer is `nr * nc` contiguous doubles,
//! column `c` is the contiguous slice `x[c*nr .. (c+1)*nr]`. A stack of `nl`
//! layers is `nl` such planes back to back. Parallelism is one rayon task per
//! output column across all layers, so the cache-friendly inner loop slides
//! the window vertically down contiguous memory.

use crate::threading::{get_num_threads, maybe_par};
use rayon::prelude::*;

/// Out-of-bounds policy at the matrix edges.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Edge {
    /// Window truncates: outside cells are simply absent from the window.
    Shrink,
    /// Indices reflect around the boundary (half-sample symmetric: -1 -> 0).
    Reflect,
    /// Indices clamp to the nearest edge cell.
    Nearest,
    /// Outside cells take this fixed value.
    Constant(f64),
}

#[derive(Clone, Copy, Debug)]
pub struct Dims {
    pub nr: usize,
    pub nc: usize,
}

/// Half-window sizes: `(wr - 1) / 2`, `(wc - 1) / 2`.
#[derive(Clone, Copy, Debug)]
pub struct Win {
    pub hr: usize,
    pub hc: usize,
}

/// Closed form of repeated half-sample symmetric reflection (period 2n).
#[inline]
fn reflect_idx(i: isize, n: usize) -> usize {
    let p = 2 * n as isize;
    let mut j = i % p;
    if j < 0 {
        j += p;
    }
    if j >= n as isize {
        j = p - 1 - j;
    }
    j as usize
}

/// Map a possibly out-of-range index to a source index. `None` means the cell
/// is outside the grid: skipped for `Shrink`, substituted by the caller for
/// `Constant`.
#[inline]
pub fn map_idx(i: isize, n: usize, edge: Edge) -> Option<usize> {
    if i >= 0 && (i as usize) < n {
        return Some(i as usize);
    }
    match edge {
        Edge::Shrink | Edge::Constant(_) => None,
        Edge::Nearest => Some(i.clamp(0, n as isize - 1) as usize),
        Edge::Reflect => Some(reflect_idx(i, n)),
    }
}

/// Running sums over the values a window actually contains.
///
/// In NA-aware mode NaNs are never pushed, so `sum` stays finite and `n`
/// counts only valid cells. In the fast path NaNs are pushed and poison
/// `sum`, which is how "any NA in window -> NA out" falls out for free:
/// reducers check `sum.is_nan()` where NaN would not propagate naturally.
#[derive(Clone, Copy, Default, Debug)]
pub struct Moments {
    pub sum: f64,
    pub sumsq: f64,
    pub n: u32,
}

impl Moments {
    #[inline]
    pub fn push(&mut self, v: f64) {
        self.sum += v;
        self.sumsq += v * v;
        self.n += 1;
    }

    #[inline]
    pub fn mean(&self) -> f64 {
        self.sum / self.n as f64
    }

    /// Population variance. The negative clamp uses a comparison rather than
    /// `f64::max` so NaN passes through (`max` would swallow it).
    #[inline]
    pub fn var_pop(&self) -> f64 {
        let n = self.n as f64;
        let v = (self.sumsq - self.sum * self.sum / n) / n;
        if v < 0.0 {
            0.0
        } else {
            v
        }
    }

    /// Sample variance (n - 1 denominator).
    #[inline]
    pub fn var_samp(&self) -> f64 {
        let n = self.n as f64;
        let v = (self.sumsq - self.sum * self.sum / n) / (n - 1.0);
        if v < 0.0 {
            0.0
        } else {
            v
        }
    }
}

#[inline]
fn push_val<const NA_AWARE: bool>(v: f64, scratch: &mut Vec<f64>, m: &mut Moments) {
    if NA_AWARE && v.is_nan() {
        return;
    }
    scratch.push(v);
    m.push(v);
}

/// Gather one window the slow, edge-checked way (border cells and any layout
/// the fast loop cannot prove in-bounds).
fn gather_border<const NA_AWARE: bool>(
    x: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    r: usize,
    c: usize,
    scratch: &mut Vec<f64>,
) -> Moments {
    scratch.clear();
    let mut m = Moments::default();
    for dc in -(w.hc as isize)..=(w.hc as isize) {
        match map_idx(c as isize + dc, d.nc, edge) {
            Some(cc) => {
                let base = cc * d.nr;
                for dr in -(w.hr as isize)..=(w.hr as isize) {
                    match map_idx(r as isize + dr, d.nr, edge) {
                        Some(rr) => push_val::<NA_AWARE>(x[base + rr], scratch, &mut m),
                        None => {
                            if let Edge::Constant(v) = edge {
                                push_val::<NA_AWARE>(v, scratch, &mut m);
                            }
                        }
                    }
                }
            }
            None => {
                // Whole window column lies outside the grid.
                if let Edge::Constant(v) = edge {
                    for _ in 0..(2 * w.hr + 1) {
                        push_val::<NA_AWARE>(v, scratch, &mut m);
                    }
                }
            }
        }
    }
    m
}

/// Process one output column. Interior rows of interior columns take the
/// bounds-check-free hot loop; everything else goes through `gather_border`.
fn scan_column<F, const NA_AWARE: bool>(
    x: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    f: &F,
    c: usize,
    out_col: &mut [f64],
    scratch: &mut Vec<f64>,
) where
    F: Fn(f64, &mut Vec<f64>, Moments) -> f64,
{
    let nr = d.nr;
    let interior_col = c >= w.hc && c + w.hc < d.nc;
    let r0 = w.hr;
    let r1 = nr.saturating_sub(w.hr);
    if interior_col && r0 < r1 {
        for r in 0..r0 {
            let m = gather_border::<NA_AWARE>(x, d, w, edge, r, c, scratch);
            out_col[r] = f(x[c * nr + r], scratch, m);
        }
        for r in r0..r1 {
            scratch.clear();
            let mut m = Moments::default();
            for wc in (c - w.hc)..=(c + w.hc) {
                let base = wc * nr + r - w.hr;
                for &v in &x[base..base + 2 * w.hr + 1] {
                    push_val::<NA_AWARE>(v, scratch, &mut m);
                }
            }
            out_col[r] = f(x[c * nr + r], scratch, m);
        }
        for r in r1..nr {
            let m = gather_border::<NA_AWARE>(x, d, w, edge, r, c, scratch);
            out_col[r] = f(x[c * nr + r], scratch, m);
        }
    } else {
        for r in 0..nr {
            let m = gather_border::<NA_AWARE>(x, d, w, edge, r, c, scratch);
            out_col[r] = f(x[c * nr + r], scratch, m);
        }
    }
}

/// Run a window reducer over a stack of layers. `x` and `out` hold
/// `k * nr * nc` doubles for some integer number of layers `k`.
///
/// `f(center, window_values, moments) -> output` is called once per cell.
/// `window_values` is a reused scratch buffer holding the gathered window
/// (NaNs excluded when `NA_AWARE`); `moments` are accumulated over the same
/// values during the gather.
pub fn run_scan<F, const NA_AWARE: bool>(
    x: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    f: &F,
    out: &mut [f64],
) where
    F: Fn(f64, &mut Vec<f64>, Moments) -> f64 + Sync,
{
    let plane = d.nr * d.nc;
    debug_assert!(plane > 0 && x.len() == out.len() && x.len() % plane == 0);
    let cap = (2 * w.hr + 1) * (2 * w.hc + 1);
    let column = |g: usize, out_col: &mut [f64], scratch: &mut Vec<f64>| {
        let (l, c) = (g / d.nc, g % d.nc);
        let layer = &x[l * plane..(l + 1) * plane];
        scan_column::<F, NA_AWARE>(layer, d, w, edge, f, c, out_col, scratch);
    };
    if get_num_threads() <= 1 {
        let mut scratch = Vec::with_capacity(cap);
        for (g, out_col) in out.chunks_mut(d.nr).enumerate() {
            column(g, out_col, &mut scratch);
        }
    } else {
        maybe_par(|| {
            out.par_chunks_mut(d.nr).enumerate().for_each_init(
                || Vec::with_capacity(cap),
                |scratch, (g, out_col)| column(g, out_col, scratch),
            );
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const NA: f64 = f64::NAN;

    fn mean_reducer(_x: f64, _v: &mut Vec<f64>, m: Moments) -> f64 {
        if m.n == 0 {
            NA
        } else {
            m.mean()
        }
    }

    fn approx(a: f64, b: f64) -> bool {
        (a.is_nan() && b.is_nan()) || (a - b).abs() < 1e-12
    }

    #[test]
    fn map_idx_policies() {
        assert_eq!(map_idx(2, 5, Edge::Shrink), Some(2));
        assert_eq!(map_idx(-1, 5, Edge::Shrink), None);
        assert_eq!(map_idx(5, 5, Edge::Constant(0.0)), None);
        assert_eq!(map_idx(-1, 5, Edge::Nearest), Some(0));
        assert_eq!(map_idx(7, 5, Edge::Nearest), Some(4));
        assert_eq!(map_idx(-1, 5, Edge::Reflect), Some(0));
        assert_eq!(map_idx(-2, 5, Edge::Reflect), Some(1));
        assert_eq!(map_idx(5, 5, Edge::Reflect), Some(4));
        assert_eq!(map_idx(6, 5, Edge::Reflect), Some(3));
        // Multi-wrap reflection for windows wider than the grid.
        assert_eq!(map_idx(-6, 5, Edge::Reflect), Some(4));
        assert_eq!(map_idx(10, 5, Edge::Reflect), Some(0));
    }

    #[test]
    fn shrink_mean_matches_hand_computed() {
        // 3x3 matrix, column-major: columns are [1,2,3], [4,5,6], [7,8,9].
        let x = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let d = Dims { nr: 3, nc: 3 };
        let w = Win { hr: 1, hc: 1 };
        let mut out = [0.0; 9];
        run_scan::<_, true>(&x, d, w, Edge::Shrink, &mean_reducer, &mut out);
        // Center = full 3x3 mean = 5; corner (0,0) = mean(1,2,4,5) = 3.
        assert!(approx(out[4], 5.0));
        assert!(approx(out[0], 3.0));
        // Edge midpoint (1,0) = mean(1,2,3,4,5,6) = 3.5.
        assert!(approx(out[1], 3.5));
    }

    #[test]
    fn constant_edge_pads_value() {
        let x = [1., 1., 1., 1.];
        let d = Dims { nr: 2, nc: 2 };
        let w = Win { hr: 1, hc: 1 };
        let mut out = [0.0; 4];
        run_scan::<_, true>(&x, d, w, Edge::Constant(10.0), &mean_reducer, &mut out);
        // Every window: 4 ones + 5 tens = 54 / 9 = 6.
        for v in out {
            assert!(approx(v, 6.0));
        }
    }

    #[test]
    fn na_aware_skips_and_fast_path_poisons() {
        let x = [1., NA, 3., 4., 5., 6., 7., 8., 9.];
        let d = Dims { nr: 3, nc: 3 };
        let w = Win { hr: 1, hc: 1 };
        let mut out = [0.0; 9];
        run_scan::<_, true>(&x, d, w, Edge::Shrink, &mean_reducer, &mut out);
        // Center: mean of the 8 valid = (1+3+4+5+6+7+8+9)/8 = 5.375.
        assert!(approx(out[4], 5.375));
        run_scan::<_, false>(&x, d, w, Edge::Shrink, &mean_reducer, &mut out);
        assert!(out[4].is_nan());
        assert!(!out[7].is_nan()); // window away from the NA is untouched
    }

    #[test]
    fn window_larger_than_grid_falls_back_to_border_path() {
        let x = [1., 2., 3., 4.];
        let d = Dims { nr: 2, nc: 2 };
        let w = Win { hr: 2, hc: 2 };
        let mut out = [0.0; 4];
        run_scan::<_, true>(&x, d, w, Edge::Shrink, &mean_reducer, &mut out);
        for v in out {
            assert!(approx(v, 2.5));
        }
        // Reflection repeats cells unevenly once the window exceeds the grid:
        // at (0,0) rows/cols carry weights (2,3), at (1,1) weights (3,2).
        run_scan::<_, true>(&x, d, w, Edge::Reflect, &mean_reducer, &mut out);
        assert!(approx(out[0], 70.0 / 25.0));
        assert!(approx(out[3], 55.0 / 25.0));
    }

    #[test]
    fn multi_layer_stack_processes_each_plane() {
        let x = [1., 1., 1., 1., 5., 5., 5., 5.];
        let d = Dims { nr: 2, nc: 2 };
        let w = Win { hr: 1, hc: 1 };
        let mut out = [0.0; 8];
        run_scan::<_, true>(&x, d, w, Edge::Shrink, &mean_reducer, &mut out);
        assert!(out[..4].iter().all(|&v| approx(v, 1.0)));
        assert!(out[4..].iter().all(|&v| approx(v, 5.0)));
    }
}
