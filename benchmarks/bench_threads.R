# Thread scaling for the two engine paths.
# Run from the package root with the package installed (release build):
#   Rscript benchmarks/bench_threads.R

library(rustyfilters)

n <- 4000L
set.seed(42)
m <- matrix(rnorm(n * n), n, n)

time_med <- function(f, reps = 3L) {
  t <- vapply(seq_len(reps), function(i) {
    system.time(f())[["elapsed"]]
  }, numeric(1))
  stats::median(t)
}

max_t <- rf_get_threads()
threads <- unique(pmin(c(1L, 2L, 4L, 8L, 16L), max_t))

cat(sprintf("matrix %d x %d\n\n", n, n))
base <- NULL
for (nt in threads) {
  rf_set_threads(nt)
  t_mean <- time_med(function() rf_mean(m, 7L, na_policy = "propagate"))
  t_med <- time_med(function() rf_median(m, 7L, na_policy = "propagate"))
  if (is.null(base)) base <- c(t_mean, t_med)
  cat(sprintf(
    "%2d threads | mean %.3fs (%.1fx) | median %.3fs (%.1fx)\n",
    nt, t_mean, base[1] / t_mean, t_med, base[2] / t_med
  ))
}
rf_set_threads(max_t)
