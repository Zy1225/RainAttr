# Permutation-Based Procedure for Rainfall Enhancement Trial Data

Implements the permutation-based procedure used in
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
for inference on attribution and SATE in rainfall enhancement trial data
under randomized assignment of ionizer operation. It allows for flexible
permutation schemes, including permutations between ionizers, between
days, or between gauge-day combinations.

## Usage

``` r
permutation_ionizer(
  B_permutation,
  permute_between_ionizer,
  permute_all_ionizers_between_day,
  permute_between_gaugeday,
  ionizer_operation_input,
  gaugeday_downwind_input,
  data_target_column_names,
  ionizer_operation_year_column_name,
  ionizer_operation_day_column_name,
  permutation_seed,
  permutation_parallel,
  permutation_parallel_num_worker,
  data,
  downwind_lmm_formula,
  downwind_propensity_formula,
  attr_type,
  x_downwind_name,
  target_only,
  rain_col_name
)
```

## Arguments

- B_permutation:

  An integer specifying the number of permutation replicates.
  (User-configurable permutation option using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- permute_between_ionizer:

  Logical. If `TRUE`, for each day, a random permutation is performed
  among the operation statuses of all ionizers that have been deployed
  on that day. (User-configurable permutation option using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- permute_all_ionizers_between_day:

  Logical. If `TRUE`, for each year, a random permutation is performed
  among the daily operation schedules of all trial days within that
  year. (User-configurable permutation option using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- permute_between_gaugeday:

  Logical. If `TRUE`, for each year, a random permutation is performed
  among the gauge-day level operation schedule of all gauge-days within
  that year. (User-configurable permutation option using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- ionizer_operation_input:

  A data frame containing ionizer operation indicators for each day
  (row) and each ionizer (column), where 1 indicates that an ionizer is
  turned on, 0 indicates that it is off and NA indicates that it is not
  deployed yet. Additionally, this data frame must include two columns
  with names specified by `ionizer_operation_day_column_name` and
  `ionizer_operation_year_column_name`, containing the day and year for
  each row. Each day must appear only once in this data frame (no
  duplicated day entries). The ionizer columns must appear in the same
  order as specified by `data_target_column_names` and must be
  consistent with the column order in `gaugeday_downwind_input`.
  (User-supplied using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- gaugeday_downwind_input:

  A binary matrix indicating which gauge-day observations (row) are
  downwind of which ionizers (column), where 1 indicates that the gauge
  is downwind of the ionizer on that day, 0 indicates that the gauge is
  not downwind of the ionizer, and NA indicates that the ionizer has not
  been deployed yet. The row order of this matrix must match that of
  `data`. The column order must correspond to the ionizer columns in
  `ionizer_operation_input` (excluding the day and year columns) and be
  in the same order as specified by `data_target_column_names`.
  (User-supplied using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- data_target_column_names:

  A character vector specifying the column names of `data` corresponding
  to the binary target indicators used in the downwind (second stage)
  LMM fitting. The order of names in this vector must match the column
  order of the corresponding ionizers in `ionizer_operation_input`
  (excluding the day and year columns) and in `gaugeday_downwind_input`.
  (User-configurable permutation option using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- ionizer_operation_year_column_name:

  A character string specifying the column name of
  `ionizer_operation_input` containing the year of each day.
  (User-configurable permutation option using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- ionizer_operation_day_column_name:

  A character string specifying the column name of
  `ionizer_operation_input` containing the day of each observation. The
  same column name should also be found in `data`. (User-configurable
  permutation option using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- permutation_seed:

  An integer specifying the random seed for the permutation-based
  procedure. Reproducibility is guaranteed only if
  `permutation_parallel` is the same, since parallel execution changes
  the order of random number generation. (User-configurable permutation
  option using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- permutation_parallel:

  Logical. If `TRUE`, each permutation run is executed in parallel
  across multiple workers. If `FALSE`, they are run sequentially.
  (User-configurable permutation option using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- permutation_parallel_num_worker:

  An integer specifying the number of parallel workers to use when
  `permutation_parallel = TRUE`. (User-configurable permutation option
  using
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md))

- data:

  A data frame containing the original dataset used in
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md),
  along with an additional column containing the fitted values generated
  from the upwind (first stage) LMM. Its column names should contain
  `data_target_column_names` (binary target indicators) and
  `ionizer_operation_day_column_name` (day). The row order of this data
  frame must match that of `gaugeday_downwind_input`. (Internal argument
  set automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- downwind_lmm_formula:

  A two sided linear formula object to be used in
  [lmer](https://rdrr.io/pkg/lme4/man/lmer.html), describing both the
  fixed-effects and random intercept part of the downwind (second stage)
  LMM. (Internal argument set automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- downwind_propensity_formula:

  A two sided linear formula object to be used in
  [`glm`](https://rdrr.io/r/stats/glm.html) with `family = "binomial"`,
  for fitting a propensity score model to the treatment indicators of
  downwind (second stage) observations. (Internal argument set
  automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- attr_type:

  A character string specifying the type of attribution estimates. Must
  be one of `"ChambersEtAl"`, `"ChambersEtAl_No_Winsorize"`,
  `"ThoEtAl"`, or `"No"`. See
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
  for more information. (Internal argument set automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- x_downwind_name:

  A character vector containing variable names from the right hand side
  of `downwind_lmm_formula`, for those variables that are not related to
  ionizers (treatment). The intercept is always included and does not
  need to be specified. (Internal argument set automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- target_only:

  Logical. If `TRUE` the attribution estimates are computed based on
  only target observations. If `FALSE` the attribution estimates are
  computed based on both treatment and control observations. (Internal
  argument set automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- rain_col_name:

  A character string specifying the column name of the raw scale
  rainfall in `ori_data`. (Internal argument set automatically when
  using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

## Value

A list with two components:

- hatattr:

  Matrix of permutation replicates of attribution estimates.

- hatsate:

  Matrix of permutation replicates of SATE estimates.

## Details

This function implements the permutation-based procedure used in
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).
It is intended for internal use only. Users should not call this
function directly. Instead, permutation-based inference should be
performed by calling
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
with `permutation = TRUE`.

To perform permutation-based procedure, additional information need to
be supplied through the following arguments of
[`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md):

- `ionizer_operation_input`:

  Day(group)-level ionizers operation schedule during the rainfall
  enhancement trial.

- gaugeday_downwind_input:

  Gauge-day(unit within group)-level information on relative orientation
  of gauges from ionizers each day.

**Permutation Steps**  
The permutation-based procedure considers to randomly permute ionizers'
operation statuses via:

1.  If `permute_between_ionizer = TRUE`, for each row of the day-level
    `ionizer_operation_input` matrix, a random permutation is performed
    among the binary indicators in the row that correspond to ionizers
    that have already been deployed during the year of the row. This is
    equivalent to randomly permuting the operation statuses of deployed
    ionizers for each day.

2.  If `permute_all_ionizers_between_day = TRUE`, for each year, a
    random permutation is performed among all rows belonging to the year
    in the day-level `ionizer_operation_input` matrix. This is
    equivalent to randomly permuting the daily operation schedules of
    all trial days belonging to each year.

3.  The permuted day-level `ionizer_operation_input` are then expanded
    into a gauge-day level binary ionizers' operation indicator matrix,
    to match the gauge-day level `data`.

4.  If `permute_between_gaugeday = TRUE`, for each year, a random
    permutation is performed among all rows belonging to the year in the
    gauge-day level binary ionizers' operation indicator matrix. This is
    equivalent to randomly permuting the gauge-day level operation
    schedules of all gauge-days belonging to each year.

The above three optional permutation steps are performed in sequence,
resulting in a final permuted gauge-day level binary ionizers' operation
indicator matrix. An elementwise multiplication is carried out between
this matrix and `gaugeday_downwind_input`, and the results are used to
replace the original columns (`data_target_column_names`) in `data` that
contain the binary target indicators used in the downwind (second stage)
LMM fitting, where NAs (for cases where the ionizer has not been
deployed yet) are replaced by zeros. Therefore, it is important to
ensure that

- The row order of `gaugeday_downwind_input` and `data` is consistent.

- The column orders of `gaugeday_downwind_input` and
  `ionizer_operation_input` is consistent, and matches the order
  specified by `data_target_column_names`.

Based on these permuted binary target indicators, the original binary
indicators \\I\_{ij}\\ for exposure to ionizers (treatment) are also
updated to be their permuted counterparts \\I\_{ij}^\*\\, where
\\I\_{ij}^\*\\ only equals zero if all permuted binary target indicators
in its corresponding row equal zero.

**Permutation Distribution of Attribution and SATE**  
The same procedure described in
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
is then repeated on the permuted dataset, where the columns containing
the binary target indicators and the binary indicators \\I\_{ij}\\ for
exposure to ionizers (treatment) are replaced according to the
permutation. This includes the fitting of the downwind (second stage)
LMM, downwind (second stage) target-only LMM, downwind (second stage)
control-only LMM to the subset of observations from `data` that are
downwind (second stage) and with positive rainfall, along with the
fitting of downwind propensity score model to the subset of observations
from `data` that are downwind (second stage) with the response being the
permuted indicator \\I\_{ij}^\*\\ for exposure to ionizers (treatment).
Finally, two attribution estimates and the SATE estimates are computed
based on the estimation results of these models fitted to the permuted
dataset, where \\z\_{ij}\\ (ionizer related covariate vector constructed
from the binary target indicators) and \\I\_{ij}\\ are replaced by their
permuted counterparts \\z\_{ij}^\*\\ and \\I\_{ij}^\*\\, respectively.

By repeatedly permuting ionizers' operation schedules, fitting models
and computing attribution and SATE estimates for `B_permutation` number
of times, this function returns the permutation distributions of

- Two attribution estimates

- Five SATE estimates

**Permutation-Based P-Value and Permutation-Based Plots**  
The permutation distributions of attribution and SATE produced by this
function are further used in
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
to:

- Compute permutation-based p-value as the proportion of permuted
  estimates that are greater than or equal to the original estimate,
  i.e., \$\$\frac{1}{B} \sum\_{b=1}^{B} 1\_{ \\ \hat{\theta}^\*\_{b}
  \geq \hat{\theta} \\ } ,\$\$ where \\\hat{\theta}^\*\_{b}\\ is the
  estimate from the \\b\\-th permuted dataset and \\\hat{\theta}\\ is
  the corresponding estimate from the original dataset.

- Plot the kernel density estimate with a solid vertical line for the
  original estimate \\\hat{\theta}\\.

## See also

[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
for the main function,
[`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md)
for specifying permutation options
