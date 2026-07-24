# Permutation Options

This function generates a list of settings controlling how
permutation-based procedures are executed within
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).

## Usage

``` r
permutation_opt(
  B_permutation = 10000,
  permute_between_ionizer = T,
  permute_all_ionizers_between_day = F,
  permute_between_gaugeday = T,
  ionizer_operation_input = ionizer_operation,
  gaugeday_downwind_input = gaugeday_downwind,
  data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03",
    "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08",
    "Target.H.09", "Target.H.10"),
  ionizer_operation_year_column_name = "Year",
  ionizer_operation_day_column_name = "TrialDay",
  permutation_seed = NULL,
  permutation_parallel = F,
  permutation_parallel_num_worker = parallel::detectCores() - 1
)
```

## Arguments

- B_permutation:

  An integer specifying the number of permutation replicates. Default is
  10000.

- permute_between_ionizer:

  Logical. If `TRUE`, for each day, a random permutation is performed
  among the operation statuses of all ionizers that have been deployed
  on that day. Default is `TRUE`.

- permute_all_ionizers_between_day:

  Logical. If `TRUE`, for each year, a random permutation is performed
  among the daily operation schedules of all trial days within that
  year. Default is `FALSE`.

- permute_between_gaugeday:

  Logical. If `TRUE`, for each year, a random permutation is performed
  among the gauge-day level operation schedule of all gauge-days within
  that year. Default is `TRUE`.

- ionizer_operation_input:

  A data frame containing ionizer operation indicators for each day
  (row) and each ionizer (column), where 1 indicates that an ionizer is
  turned on, 0 indicates that it is off, and NA indicates that it is not
  deployed yet. Additionally, this data frame must include two columns
  with names specified by `ionizer_operation_day_column_name` and
  `ionizer_operation_year_column_name`, containing the day and year for
  each row. Each day must appear only once in this data frame (no
  duplicated day entries). The ionizer columns must appear in the same
  order as specified by `data_target_column_names` and must be
  consistent with the column order in `gaugeday_downwind_input`. Default
  is `ionizer_operation`.

- gaugeday_downwind_input:

  A binary matrix indicating which gauge-day observations (row) are
  downwind of which ionizers (column), where 1 indicates that the gauge
  is downwind of the ionizer on that day, 0 indicates that the gauge is
  not downwind of the ionizer, and NA indicates that the ionizer has not
  been deployed yet. The row order of this matrix must match that of the
  original dataset supplied to
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).
  The column order must correspond to the ionizer columns in
  `ionizer_operation_input` (excluding the day and year columns) and be
  in the same order as specified by `data_target_column_names`. Default
  is `gaugeday_downwind`.

- data_target_column_names:

  A character vector specifying the column names of the original dataset
  supplied to
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md),
  corresponding to the binary target indicators used in the downwind
  (second stage) LMM fitting. The order of names in this vector must
  match the column order of the corresponding ionizers in
  `ionizer_operation_input` (excluding the day and year columns) and in
  `gaugeday_downwind_input`. Default is:
  `c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10")`.

- ionizer_operation_year_column_name:

  A character string specifying the column name of
  `ionizer_operation_input` containing the year of each day. Default is
  `'Year'`.

- ionizer_operation_day_column_name:

  A character string specifying the column name of
  `ionizer_operation_input` containing the day of each observation. The
  same column name should also be found in the original dataset supplied
  to
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).
  Default is `'TrialDay'`.

- permutation_seed:

  An integer specifying the random seed for the permutation-based
  procedure. Reproducibility is guaranteed only if
  `permutation_parallel` is the same, since parallel execution changes
  the order of random number generation. Default is `NULL`, meaning no
  seed is set internally and users should call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand to
  ensure reproducibility.

- permutation_parallel:

  Logical. If `TRUE`, each permutation run is executed in parallel
  across multiple workers. If `FALSE`, they are run sequentially.
  Default is `FALSE`.

- permutation_parallel_num_worker:

  An integer specifying the number of parallel workers to use when
  `permutation_parallel = TRUE`. Default is
  `parallel::detectCores() - 1`.

## Value

A list containing all permutation options, suitable for passing to
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).

## Details

This function is used to configure and store all settings needed for
performing permutation-based analyses within
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).
It is important to ensure that

- The row order of `gaugeday_downwind_input` and the original dataset
  supplied to
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
  is consistent.

- The column orders of `gaugeday_downwind_input` and
  `ionizer_operation_input` is consistent, and matches the order
  specified by `data_target_column_names`.

## See also

[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
for the main function,
[`permutation_ionizer`](https://zy1225.github.io/RainAttr/reference/permutation_ionizer.md)
for more details on the permutation-based procedure

## Examples

``` r
#Create default permutation options
# These are the same permutation settings used in Chambers et al. (2022a)
# Nudging a Pseudo-Science Towards a Science—The Role of Statistics
# in a Rainfall Enhancement Trial in Oman. \emph{International Statistical Review}, 90: 346–373,
# as well as Chambers et al. (2022b) Weighting, Informativeness and Causal Inference,
# with an Application to Rainfall Enhancement. \emph{Journal of the Royal Statistical Society Series A: Statistics in Society}, 185: 1584–1612.
# Specifically: permute_between_ionizer = TRUE, permute_all_ionizers_between_day = FALSE, and permute_between_gaugeday = TRUE
perm_options = permutation_opt()
str(perm_options)
#> List of 12
#>  $ B_permutation                     : num 10000
#>  $ permute_between_ionizer           : logi TRUE
#>  $ permute_all_ionizers_between_day  : logi FALSE
#>  $ permute_between_gaugeday          : logi TRUE
#>  $ ionizer_operation_input           :'data.frame':  740 obs. of  12 variables:
#>   ..$ TrialDay: num [1:740] 2013135 2013136 2013137 2013138 2013139 ...
#>   ..$ Year    : num [1:740] 2013 2013 2013 2013 2013 ...
#>   ..$ H1      : num [1:740] 1 1 0 0 0 1 0 1 0 0 ...
#>   ..$ H2      : num [1:740] 0 0 1 1 1 0 1 0 1 1 ...
#>   ..$ H3      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H4      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H5      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H6      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H7      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H8      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H9      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H10     : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>  $ gaugeday_downwind_input           : num [1:122259, 1:10] 0 0 0 0 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:10] "H1" "H2" "H3" "H4" ...
#>  $ data_target_column_names          : chr [1:10] "Target.H.01" "Target.H.02" "Target.H.03" "Target.H.04" ...
#>  $ ionizer_operation_year_column_name: chr "Year"
#>  $ ionizer_operation_day_column_name : chr "TrialDay"
#>  $ permutation_seed                  : NULL
#>  $ permutation_parallel              : logi FALSE
#>  $ permutation_parallel_num_worker   : num 3

#Permutation option with parallelization over (parallel::detectCores() - 1) number of workers and seed = 1 for reproducibility
perm_options_parallel = permutation_opt(
  permutation_seed = 1,
  permutation_parallel = TRUE,
  permutation_parallel_num_worker = parallel::detectCores() - 1
)
str(perm_options_parallel)
#> List of 12
#>  $ B_permutation                     : num 10000
#>  $ permute_between_ionizer           : logi TRUE
#>  $ permute_all_ionizers_between_day  : logi FALSE
#>  $ permute_between_gaugeday          : logi TRUE
#>  $ ionizer_operation_input           :'data.frame':  740 obs. of  12 variables:
#>   ..$ TrialDay: num [1:740] 2013135 2013136 2013137 2013138 2013139 ...
#>   ..$ Year    : num [1:740] 2013 2013 2013 2013 2013 ...
#>   ..$ H1      : num [1:740] 1 1 0 0 0 1 0 1 0 0 ...
#>   ..$ H2      : num [1:740] 0 0 1 1 1 0 1 0 1 1 ...
#>   ..$ H3      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H4      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H5      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H6      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H7      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H8      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H9      : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>   ..$ H10     : num [1:740] NA NA NA NA NA NA NA NA NA NA ...
#>  $ gaugeday_downwind_input           : num [1:122259, 1:10] 0 0 0 0 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:10] "H1" "H2" "H3" "H4" ...
#>  $ data_target_column_names          : chr [1:10] "Target.H.01" "Target.H.02" "Target.H.03" "Target.H.04" ...
#>  $ ionizer_operation_year_column_name: chr "Year"
#>  $ ionizer_operation_day_column_name : chr "TrialDay"
#>  $ permutation_seed                  : num 1
#>  $ permutation_parallel              : logi TRUE
#>  $ permutation_parallel_num_worker   : num 3
```
