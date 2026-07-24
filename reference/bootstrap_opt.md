# Bootstrap Options

This function generates a list of settings controlling how two-level
bootstrap procedures are executed within
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).

## Usage

``` r
bootstrap_opt(
  B_bootstrap = 10000,
  bootstrap_type = "PREB1",
  bootstrap_zero = T,
  positive_prob_threshold = NULL,
  discretize_rain = T,
  winsorize_individual_rain = T,
  individual_rain_interval = c(100, 175),
  winsorize_total_rain = T,
  total_rain_interval = c(6000, 60000),
  bootstrap_seed = NULL,
  bootstrap_parallel = F,
  bootstrap_parallel_num_worker = parallel::detectCores() - 1,
  CI_level = 0.95
)
```

## Arguments

- B_bootstrap:

  An integer specifying the number of bootstrap replicates. Default is
  10000.

- bootstrap_type:

  A character string specifying the type of bootstrap. Must be one of
  `"REB0"`, `"REB1"`, `"REB2"`, `"PREB0"`, `"PREB1"`, `"PREB2"`, or
  `"MREB1"`. See
  [`bootstrap_downwind`](https://zy1225.github.io/RainAttr/reference/bootstrap_downwind.md)
  for their differences. Default is "PREB1".

- bootstrap_zero:

  Logical. If `TRUE`, the optional first-level bootstrap is performed to
  generate bootstrap samples of binary rainfall event indicators.
  Default is `TRUE`.

- positive_prob_threshold:

  An optional numeric value between 0 and 1 specifying the probability
  threshold for generating bootstrap samples of binary rainfall event
  indicators. Probabilities below this threshold are set to zero.
  Default is `NULL`.

- discretize_rain:

  Logical. If `TRUE`, rainfall values are discretized in bootstrap
  resamples. Default is `TRUE`.

- winsorize_individual_rain:

  Logical. If `TRUE`, individual rainfall values in bootstrap samples
  that exceed the upper bound specified by `individual_rain_interval`
  are replaced with random draws from a uniform distribution over
  `[individual_rain_interval[1], individual_rain_interval[2]]`. Default
  is `TRUE`.

- individual_rain_interval:

  Numeric vector of length 2 specifying the lower and upper bounds for
  adjusting bootstrapped individual rainfall values that are too large
  when `winsorize_individual_rain = TRUE`. Default is `c(100,175)`.

- winsorize_total_rain:

  Logical. If `TRUE`, all individual rainfall values in each bootstrap
  sample are proportionally rescaled so that the total equals a random
  number drawn uniformly from
  `[total_rain_interval[1], total_rain_interval[2]]` whenever the total
  bootstrapped rainfall falls outside this interval. Default is `TRUE`.

- total_rain_interval:

  Numeric vector of length 2 specifying the lower and upper bounds for
  adjusting the total of bootstrapped rainfall values when
  `winsorize_total_rain = TRUE`. Default is `c(6000,60000)`.

- bootstrap_seed:

  An integer specifying the random seed for the bootstrap procedure.
  Reproducibility is guaranteed only if `bootstrap_parallel` is the
  same, since parallel execution changes the order of random number
  generation. Default is `NULL`, meaning no seed is set internally and
  users should call [`set.seed()`](https://rdrr.io/r/base/Random.html)
  beforehand to ensure reproducibility.

- bootstrap_parallel:

  Logical. If `TRUE`, each bootstrap run is executed in parallel across
  multiple workers. If `FALSE`, they are run sequentially. Default is
  `FALSE`.

- bootstrap_parallel_num_worker:

  An integer specifying the number of parallel workers to use when
  `bootstrap_parallel = TRUE`. Default is `parallel::detectCores() - 1`.

- CI_level:

  A numeric value between 0 and 1 specifying the confidence level of the
  bootstrap percentile confidence intervals. Default is 0.95.

## Value

A list containing all bootstrap options, suitable for passing to
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).

## Details

This function is used to configure and store all settings needed for
performing two-level bootstrap analyses within
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).

## See also

[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
for the main function,
[`bootstrap_downwind`](https://zy1225.github.io/RainAttr/reference/bootstrap_downwind.md)
for more details on the two-level bootstrap procedure

## Examples

``` r
#Create default bootstrap options to account for highly unbalanced clustered data
# Specifically: bootstrap_type = 'PREB1' as proposed by Tho et al. (2025) Adjusted Random Effect Block Bootstraps for Highly Unbalanced Clustered Data. arXiv:2510.07770.
boot_options = bootstrap_opt()
str(boot_options)
#> List of 13
#>  $ B_bootstrap                  : num 10000
#>  $ bootstrap_type               : chr "PREB1"
#>  $ bootstrap_zero               : logi TRUE
#>  $ positive_prob_threshold      : NULL
#>  $ discretize_rain              : logi TRUE
#>  $ winsorize_individual_rain    : logi TRUE
#>  $ individual_rain_interval     : num [1:2] 100 175
#>  $ winsorize_total_rain         : logi TRUE
#>  $ total_rain_interval          : num [1:2] 6000 60000
#>  $ bootstrap_seed               : NULL
#>  $ bootstrap_parallel           : logi FALSE
#>  $ bootstrap_parallel_num_worker: num 3
#>  $ CI_level                     : num 0.95

#Bootstrap option with parallelization over
#(parallel::detectCores() - 1) number of workers
#and seed = 1 for reproducibility
boot_options_parallel = bootstrap_opt(
  bootstrap_seed = 1,
  bootstrap_parallel = TRUE,
  bootstrap_parallel_num_worker = parallel::detectCores() - 1
)
str(boot_options_parallel)
#> List of 13
#>  $ B_bootstrap                  : num 10000
#>  $ bootstrap_type               : chr "PREB1"
#>  $ bootstrap_zero               : logi TRUE
#>  $ positive_prob_threshold      : NULL
#>  $ discretize_rain              : logi TRUE
#>  $ winsorize_individual_rain    : logi TRUE
#>  $ individual_rain_interval     : num [1:2] 100 175
#>  $ winsorize_total_rain         : logi TRUE
#>  $ total_rain_interval          : num [1:2] 6000 60000
#>  $ bootstrap_seed               : num 1
#>  $ bootstrap_parallel           : logi TRUE
#>  $ bootstrap_parallel_num_worker: num 3
#>  $ CI_level                     : num 0.95
```
