# Benchmark rustyfilters against terra::focal on a large matrix.
# Run from the package root with the package installed (release build):
#   Rscript benchmarks/bench_filters.R

library(rustyfilters)

n <- 4000L
set.seed(42)
m <- matrix(rnorm(n * n), n, n)

time_med <- function(expr, reps = 3L) {
  t <- vapply(seq_len(reps), function(i) {
    system.time(force(expr))[["elapsed"]]
  }, numeric(1))
  stats::median(t)
}

cat(sprintf("matrix %d x %d, %d threads\n\n", n, n, rf_get_threads()))

for (w in c(3L, 7L, 15L)) {
  t_mean <- time_med(rf_mean(m, w, na_policy = "propagate"))
  t_sd <- time_med(rf_focal(m, w, "sd", na_policy = "propagate"))
  t_med <- time_med(rf_median(m, w, na_policy = "propagate"))
  t_lee <- time_med(rf_lee(m, w))
  cat(sprintf(
    "window %2d | mean %.3fs | sd %.3fs | median %.3fs | lee %.3fs\n",
    w, t_mean, t_sd, t_med, t_lee
  ))
}

if (requireNamespace("terra", quietly = TRUE)) {
  r <- terra::rast(m)
  cat("\nterra::focal comparison (window 7):\n")
  t_ours <- time_med(rf_mean(m, 7L, na_policy = "propagate"))
  t_terra <- time_med(terra::focal(r, 7, "mean", na.rm = FALSE))
  cat(sprintf(
    "mean   | rustyfilters %.3fs | terra %.3fs | %.1fx\n",
    t_ours, t_terra, t_terra / t_ours
  ))
  t_ours <- time_med(rf_median(m, 7L, na_policy = "propagate"))
  t_terra <- time_med(terra::focal(r, 7, "median", na.rm = FALSE))
  cat(sprintf(
    "median | rustyfilters %.3fs | terra %.3fs | %.1fx\n",
    t_ours, t_terra, t_terra / t_ours
  ))
}
