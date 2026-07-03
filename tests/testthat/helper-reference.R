# Naive pure-R moving-window reference: the oracle every Rust filter is
# tested against. Slow double loop, but obviously correct.

# Map a possibly out-of-range 1-based index onto 1..n, or NA when the cell is
# outside the grid (shrink/constant).
ref_map_idx <- function(i, n, edge) {
  if (i >= 1 && i <= n) {
    return(i)
  }
  switch(edge,
    shrink = ,
    constant = NA_integer_,
    nearest = min(max(i, 1L), n),
    reflect = {
      j <- (i - 1L) %% (2L * n)
      if (j < 0L) j <- j + 2L * n
      if (j >= n) j <- 2L * n - 1L - j
      j + 1L
    }
  )
}

# Gather the window values around (r, c) under an edge policy. Outside cells
# are dropped (shrink) or replaced by edge_value (constant).
ref_window <- function(x, r, c, w, edge, edge_value) {
  hr <- (w[1] - 1) / 2
  hc <- (w[2] - 1) / 2
  vals <- numeric(0)
  for (dc in -hc:hc) {
    cc <- ref_map_idx(c + dc, ncol(x), edge)
    for (dr in -hr:hr) {
      rr <- ref_map_idx(r + dr, nrow(x), edge)
      if (is.na(cc) || is.na(rr)) {
        if (edge == "constant") vals <- c(vals, edge_value)
      } else {
        vals <- c(vals, x[rr, cc])
      }
    }
  }
  vals
}

ref_stat <- function(vals, stat, na_policy) {
  if (na_policy == "omit") {
    vals <- vals[!is.na(vals)]
  } else if (anyNA(vals)) {
    return(NA_real_)
  }
  if (length(vals) == 0) {
    return(NA_real_)
  }
  switch(stat,
    mean = mean(vals),
    median = stats::median(vals),
    min = min(vals),
    max = max(vals),
    range = max(vals) - min(vals),
    sd = if (length(vals) < 2) NA_real_ else stats::sd(vals),
    sum = sum(vals),
    mode = {
      u <- sort(unique(vals))
      u[which.max(vapply(u, function(z) sum(vals == z), integer(1)))]
    }
  )
}

ref_focal <- function(x, window = 3L, stat = "mean", edge = "shrink",
                      edge_value = 0, na_policy = "omit") {
  w <- rep(as.integer(window), length.out = 2L)
  out <- matrix(NA_real_, nrow(x), ncol(x), dimnames = dimnames(x))
  for (c in seq_len(ncol(x))) {
    for (r in seq_len(nrow(x))) {
      vals <- ref_window(x, r, c, w, edge, edge_value)
      out[r, c] <- ref_stat(vals, stat, na_policy)
    }
  }
  out
}

# A reproducible test matrix with optional NA holes.
ref_mat <- function(nr = 7, nc = 9, na = 0, seed = 42) {
  withr::with_seed(seed, {
    m <- matrix(round(stats::rnorm(nr * nc), 2), nr, nc)
    if (na > 0) m[sample.int(length(m), na)] <- NA_real_
    m
  })
}
