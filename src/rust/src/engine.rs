//! Core moving-window machinery shared by every filter. Pure Rust (no
//! extendr types) so it is unit-testable with plain `cargo test`.
//!
//! Data is column-major (R layout): a layer is `nr * nc` contiguous doubles,
//! column `c` is the contiguous slice `x[c*nr .. (c+1)*nr]`. A stack of `nl`
//! layers is `nl` such planes back to back. Parallelism is one rayon task per
//! output column across all layers, so the cache-friendly inner loop slides
//! the window vertically down contiguous memory.

use crate::threading::maybe_par;
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
#[allow(clippy::too_many_arguments)]
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
    maybe_par(|| {
        out.par_chunks_mut(d.nr).enumerate().for_each_init(
            || Vec::with_capacity(cap),
            |scratch, (g, out_col)| column(g, out_col, scratch),
        );
    });
}

/// Drive a per-column worker over a stack of output columns, sequentially or
/// on the rayon pool depending on the thread setting. `f(g, out_col)` gets
/// the global column index `g` (layer = g / nc, column = g % nc).
pub fn par_columns<F>(out: &mut [f64], nr: usize, f: F)
where
    F: Fn(usize, &mut [f64]) + Sync,
{
    maybe_par(|| {
        out.par_chunks_mut(nr)
            .enumerate()
            .for_each(|(g, out_col)| f(g, out_col));
    });
}

// ---------------------------------------------------------------------------
// Moments fast path: separable sliding box sums, O(1) per pixel regardless of
// window size. Used by every filter that only needs the window mean/variance
// (mean, sum, sd, Kuan, enhanced Lee, Gamma-MAP).
//
// Per rayon task: a ring buffer holds the vertical window sums of the last
// `wc` source columns; an accumulator column holds the horizontal total. Both
// slide by add/subtract, so drift is bounded by the column length and the
// band width. NaNs are never added to the sliding sums (a NaN entering a
// running sum would poison every window after it); instead the fast path
// keeps a sliding count of NaNs per window and poisons exactly the windows
// that contain one, preserving `propagate` locality.

/// Per-window vertical sums for one source column.
struct ColSums {
    sum: Vec<f64>,
    sumsq: Vec<f64>,
    cnt: Vec<u32>,
    nan: Vec<u32>,
}

impl ColSums {
    fn new(nr: usize) -> Self {
        ColSums {
            sum: vec![0.0; nr],
            sumsq: vec![0.0; nr],
            cnt: vec![0; nr],
            nan: vec![0; nr],
        }
    }
}

/// Add one cell into a running (sum, sumsq, valid count, NaN count).
#[inline]
fn acc_cell(v: f64, s: &mut f64, q: &mut f64, n: &mut u32, nn: &mut u32) {
    if v.is_nan() {
        *nn += 1;
    } else {
        *s += v;
        *q += v * v;
        *n += 1;
    }
}

/// Fill `cs` with the vertical window sums of window column `j` (possibly out
/// of grid range; edge policy decides what it holds).
fn vertical_sums<const NA_AWARE: bool>(
    x: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    j: isize,
    cs: &mut ColSums,
) {
    let nr = d.nr;
    let wr = 2 * w.hr + 1;
    let cc = match map_idx(j, d.nc, edge) {
        Some(cc) => cc,
        None => {
            // Whole column outside the grid: a constant pad or nothing.
            match edge {
                Edge::Constant(v) if !v.is_nan() => {
                    cs.sum.fill(wr as f64 * v);
                    cs.sumsq.fill(wr as f64 * v * v);
                    cs.cnt.fill(wr as u32);
                    cs.nan.fill(0);
                }
                Edge::Constant(_) => {
                    // NaN pad: all cells missing.
                    cs.sum.fill(0.0);
                    cs.sumsq.fill(0.0);
                    cs.cnt.fill(0);
                    cs.nan.fill(wr as u32);
                }
                _ => {
                    cs.sum.fill(0.0);
                    cs.sumsq.fill(0.0);
                    cs.cnt.fill(0);
                    cs.nan.fill(0);
                }
            }
            return;
        }
    };
    let col = &x[cc * nr..(cc + 1) * nr];
    let direct = |r: usize, cs: &mut ColSums| {
        let (mut s, mut q, mut n, mut nn) = (0.0, 0.0, 0u32, 0u32);
        for dr in -(w.hr as isize)..=(w.hr as isize) {
            match map_idx(r as isize + dr, nr, edge) {
                Some(rr) => acc_cell(col[rr], &mut s, &mut q, &mut n, &mut nn),
                None => {
                    if let Edge::Constant(v) = edge {
                        acc_cell(v, &mut s, &mut q, &mut n, &mut nn);
                    }
                }
            }
        }
        cs.sum[r] = s;
        cs.sumsq[r] = q;
        cs.cnt[r] = n;
        cs.nan[r] = nn;
    };
    let r1 = nr.saturating_sub(w.hr);
    if w.hr >= r1 {
        // Window taller than the column: every row is a border row.
        for r in 0..nr {
            direct(r, cs);
        }
        return;
    }
    for r in 0..w.hr {
        direct(r, cs);
    }
    // Interior rows slide: enter col[r + hr], leave col[r - hr - 1].
    let (mut s, mut q, mut n, mut nn) = (0.0, 0.0, 0u32, 0u32);
    for &v in col[0..wr].iter() {
        acc_cell(v, &mut s, &mut q, &mut n, &mut nn);
    }
    cs.sum[w.hr] = s;
    cs.sumsq[w.hr] = q;
    cs.cnt[w.hr] = n;
    cs.nan[w.hr] = nn;
    for r in (w.hr + 1)..r1 {
        let ev = col[r + w.hr];
        let lv = col[r - w.hr - 1];
        if ev.is_nan() {
            nn += 1;
        } else {
            s += ev;
            q += ev * ev;
            n += 1;
        }
        if lv.is_nan() {
            nn -= 1;
        } else {
            s -= lv;
            q -= lv * lv;
            n -= 1;
        }
        cs.sum[r] = s;
        cs.sumsq[r] = q;
        cs.cnt[r] = n;
        cs.nan[r] = nn;
    }
    for r in r1..nr {
        direct(r, cs);
    }
}

/// Turn accumulated sums into the `Moments` a reducer sees. In the fast path
/// a window containing any NaN is poisoned; in NA-aware mode NaNs were simply
/// excluded.
#[inline]
fn materialize<const NA_AWARE: bool>(s: f64, q: f64, n: u32, nn: u32) -> Moments {
    if !NA_AWARE && nn > 0 {
        Moments {
            sum: f64::NAN,
            sumsq: f64::NAN,
            n: n + nn,
        }
    } else {
        Moments {
            sum: s,
            sumsq: q,
            n,
        }
    }
}

/// Process one band of output columns `[c0, c0 + width)` of one layer.
fn moments_band<F, const NA_AWARE: bool>(
    x: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    f: &F,
    c0: usize,
    out_band: &mut [f64],
) where
    F: Fn(f64, Moments) -> f64,
{
    let nr = d.nr;
    let wc = 2 * w.hc + 1;
    let width = out_band.len() / nr;
    let mut ring: Vec<ColSums> = (0..wc).map(|_| ColSums::new(nr)).collect();
    let mut asum = vec![0.0; nr];
    let mut asumsq = vec![0.0; nr];
    let mut acnt = vec![0u32; nr];
    let mut anan = vec![0u32; nr];
    // Prime the ring with the window columns of output column c0.
    for (slot, dc) in (-(w.hc as isize)..=(w.hc as isize)).enumerate() {
        vertical_sums::<NA_AWARE>(x, d, w, edge, c0 as isize + dc, &mut ring[slot]);
        let cs = &ring[slot];
        for r in 0..nr {
            asum[r] += cs.sum[r];
            asumsq[r] += cs.sumsq[r];
            acnt[r] += cs.cnt[r];
            anan[r] += cs.nan[r];
        }
    }
    for (i, out_col) in out_band.chunks_mut(nr).enumerate() {
        let c = c0 + i;
        if i > 0 {
            // Slide: the slot holding window column (c - 1 - hc) is replaced
            // by window column (c + hc), which lands on the same ring slot.
            let slot = (i - 1) % wc;
            {
                let cs = &ring[slot];
                for r in 0..nr {
                    asum[r] -= cs.sum[r];
                    asumsq[r] -= cs.sumsq[r];
                    acnt[r] -= cs.cnt[r];
                    anan[r] -= cs.nan[r];
                }
            }
            vertical_sums::<NA_AWARE>(x, d, w, edge, (c + w.hc) as isize, &mut ring[slot]);
            let cs = &ring[slot];
            for r in 0..nr {
                asum[r] += cs.sum[r];
                asumsq[r] += cs.sumsq[r];
                acnt[r] += cs.cnt[r];
                anan[r] += cs.nan[r];
            }
        }
        let center = &x[c * nr..(c + 1) * nr];
        for r in 0..nr {
            let m = materialize::<NA_AWARE>(asum[r], asumsq[r], acnt[r], anan[r]);
            out_col[r] = f(center[r], m);
        }
    }
    debug_assert!(width >= 1);
}

/// Run a moments reducer over a stack of layers via separable sliding sums.
/// `f(center, moments) -> output` is called once per cell; cost per pixel is
/// independent of the window size.
pub fn run_moments<F, const NA_AWARE: bool>(
    x: &[f64],
    d: Dims,
    w: Win,
    edge: Edge,
    f: &F,
    out: &mut [f64],
) where
    F: Fn(f64, Moments) -> f64 + Sync,
{
    let plane = d.nr * d.nc;
    debug_assert!(plane > 0 && x.len() == out.len() && x.len() % plane == 0);
    // Fixed band width, independent of the thread count: the summation order
    // is then identical however the bands are scheduled, so results are
    // bitwise reproducible across thread settings. 64 columns keeps the ring
    // re-priming overhead small while bounding the horizontal slide's
    // floating-point drift.
    let band = 64.min(d.nc);
    for (layer_x, layer_out) in x.chunks(plane).zip(out.chunks_mut(plane)) {
        maybe_par(|| {
            layer_out
                .par_chunks_mut(band * d.nr)
                .enumerate()
                .for_each(|(b, out_band)| {
                    moments_band::<F, NA_AWARE>(layer_x, d, w, edge, f, b * band, out_band);
                });
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

    /// Deterministic pseudo-random doubles in [0, 1) with NA holes.
    fn lcg_data(n: usize, na_every: usize) -> Vec<f64> {
        let mut state: u64 = 0x2545F4914F6CDD1D;
        (0..n)
            .map(|i| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                if na_every > 0 && i % na_every == na_every - 1 {
                    NA
                } else {
                    (state >> 11) as f64 / (1u64 << 53) as f64
                }
            })
            .collect()
    }

    #[test]
    fn moments_path_matches_scan_path() {
        let d = Dims { nr: 23, nc: 17 };
        let x = lcg_data(d.nr * d.nc, 7);
        let w = Win { hr: 2, hc: 1 };
        let edges = [
            Edge::Shrink,
            Edge::Reflect,
            Edge::Nearest,
            Edge::Constant(0.5),
            Edge::Constant(NA),
        ];
        let mean_m = |_: f64, m: Moments| if m.n == 0 { NA } else { m.mean() };
        let sd_m = |_: f64, m: Moments| if m.n < 2 { NA } else { m.var_samp().sqrt() };
        let mean_s = |_: f64, _: &mut Vec<f64>, m: Moments| if m.n == 0 { NA } else { m.mean() };
        let sd_s =
            |_: f64, _: &mut Vec<f64>, m: Moments| if m.n < 2 { NA } else { m.var_samp().sqrt() };
        let mut a = vec![0.0; x.len()];
        let mut b = vec![0.0; x.len()];
        for edge in edges {
            run_moments::<_, true>(&x, d, w, edge, &mean_m, &mut a);
            run_scan::<_, true>(&x, d, w, edge, &mean_s, &mut b);
            assert!(
                a.iter().zip(&b).all(|(p, q)| approx(*p, *q)),
                "mean omit {edge:?}"
            );
            run_moments::<_, false>(&x, d, w, edge, &mean_m, &mut a);
            run_scan::<_, false>(&x, d, w, edge, &mean_s, &mut b);
            assert!(
                a.iter().zip(&b).all(|(p, q)| approx(*p, *q)),
                "mean prop {edge:?}"
            );
            run_moments::<_, true>(&x, d, w, edge, &sd_m, &mut a);
            run_scan::<_, true>(&x, d, w, edge, &sd_s, &mut b);
            assert!(
                a.iter().zip(&b).all(|(p, q)| approx(*p, *q)),
                "sd omit {edge:?}"
            );
        }
    }

    #[test]
    fn moments_path_handles_oversized_windows() {
        let d = Dims { nr: 4, nc: 3 };
        let x = lcg_data(d.nr * d.nc, 5);
        let w = Win { hr: 4, hc: 3 };
        let mean_m = |_: f64, m: Moments| if m.n == 0 { NA } else { m.mean() };
        let mean_s = |_: f64, _: &mut Vec<f64>, m: Moments| if m.n == 0 { NA } else { m.mean() };
        let mut a = vec![0.0; x.len()];
        let mut b = vec![0.0; x.len()];
        for edge in [
            Edge::Shrink,
            Edge::Reflect,
            Edge::Nearest,
            Edge::Constant(1.0),
        ] {
            run_moments::<_, true>(&x, d, w, edge, &mean_m, &mut a);
            run_scan::<_, true>(&x, d, w, edge, &mean_s, &mut b);
            assert!(a.iter().zip(&b).all(|(p, q)| approx(*p, *q)), "{edge:?}");
        }
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
