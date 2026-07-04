# Control the number of threads used by rustyfilters

All filters process windows in parallel via rayon. On load the thread
count defaults to all available cores, unless overridden by
`option(rustyfilters.threads = n)` or the `RUSTYFILTERS_NUM_THREADS`
environment variable (checked in that order), or capped at 2 when R is
running checks with `_R_CHECK_LIMIT_CORES_` set. Set to 1 for sequential
processing with zero threading overhead.

## Usage

``` r
rf_set_threads(n = 1L)

rf_get_threads()
```

## Arguments

- n:

  Integer scalar. Number of threads. Must be \>= 1.

## Value

`rf_set_threads()` invisibly returns the previous thread count;
`rf_get_threads()` returns the current thread count as an integer
scalar.

## Examples

``` r
old <- rf_set_threads(2L)
rf_get_threads()
#> [1] 2
rf_set_threads(old)
```
