# Bootstrap Procedure for Rainfall Enhancement Trial Data

Implements the two-level bootstrap procedure used in
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
for inference on attribution and sample average treatment effect in
rainfall enhancement trial data.

## Usage

``` r
bootstrap_downwind(
  B_bootstrap,
  bootstrap_type,
  bootstrap_zero,
  positive_prob_threshold = NULL,
  discretize_rain,
  winsorize_individual_rain,
  individual_rain_interval,
  winsorize_total_rain,
  total_rain_interval,
  bootstrap_seed,
  bootstrap_parallel,
  bootstrap_parallel_num_worker,
  ori_data,
  downwind,
  ori_positive,
  rain_col_name,
  downwind_target_expr,
  downwind_control_expr,
  ori_fitted_models,
  downwind_lmm_formula,
  attr_type,
  x_downwind_name,
  target_only,
  downwind_propensity_formula,
  ori_attr_est,
  ori_sate_est
)
```

## Arguments

- B_bootstrap:

  An integer specifying the number of bootstrap replicates.
  (User-configurable bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- bootstrap_type:

  A character string specifying the type of bootstrap. Must be one of
  `"REB0"`, `"REB1"`, `"REB2"`, `"PREB0"`, `"PREB1"`, `"PREB2"`, or
  `"MREB1"`. (User-configurable bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- bootstrap_zero:

  Logical. If `TRUE`, the optional first-level bootstrap is performed to
  generate bootstrap samples of binary rainfall event indicators.
  (User-configurable bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- positive_prob_threshold:

  An optional numeric value between 0 and 1 specifying the probability
  threshold for generating bootstrap samples of binary rainfall event
  indicators. Probabilities below this threshold are set to zero.
  (User-configurable bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- discretize_rain:

  Logical. If `TRUE`, rainfall values are discretized in bootstrap
  resamples. (User-configurable bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- winsorize_individual_rain:

  Logical. If `TRUE`, individual rainfall values in bootstrap samples
  that exceed the upper bound specified by `individual_rain_interval`
  are replaced with random draws from a uniform distribution over
  `[individual_rain_interval[1], individual_rain_interval[2]]`.
  (User-configurable bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- individual_rain_interval:

  Numeric vector of length 2 specifying the lower and upper bounds for
  adjusting bootstrapped individual rainfall values that are too large
  when `winsorize_individual_rain = TRUE`. (User-configurable bootstrap
  option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- winsorize_total_rain:

  Logical. If `TRUE`, all individual rainfall values in each bootstrap
  sample are proportionally rescaled so that the total equals a random
  number drawn uniformly from
  `[total_rain_interval[1], total_rain_interval[2]]` whenever the total
  bootstrapped rainfall falls outside this interval. (User-configurable
  bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- total_rain_interval:

  Numeric vector of length 2 specifying the lower and upper bounds for
  adjusting the total of bootstrapped rainfall values when
  `winsorize_total_rain = TRUE`. (User-configurable bootstrap option
  using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- bootstrap_seed:

  An integer specifying the random seed for the bootstrap procedure.
  Reproducibility is guaranteed only if `bootstrap_parallel` is the
  same, since parallel execution changes the order of random number
  generation. (User-configurable bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- bootstrap_parallel:

  Logical. If `TRUE`, each bootstrap run is executed in parallel across
  multiple workers. If `FALSE`, they are run sequentially.
  (User-configurable bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- bootstrap_parallel_num_worker:

  An integer specifying the number of parallel workers to use when
  `bootstrap_parallel = TRUE`. (User-configurable bootstrap option using
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md))

- ori_data:

  A data frame containing the original dataset used in
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md),
  along with an additional column containing the fitted values generated
  from the upwind (first stage) LMM. (Internal argument set
  automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- downwind:

  A logical vector indicating which observation in `ori_data` would be
  used in the downwind (second stage) LMM fitting. (Internal argument
  set automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- ori_positive:

  A logical vector indicating which observation in `ori_data` has
  positive rainfall. (Internal argument set automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- rain_col_name:

  A character string specifying the column name of the raw scale
  rainfall in `ori_data`. (Internal argument set automatically when
  using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- downwind_target_expr:

  A quosure (created using
  [`rlang::enquo()`](https://rlang.r-lib.org/reference/enquo.html))
  representing a logical expression used to extract the relevant subset
  of downwind (second stage) observations from `ori_data` that were
  exposed to treatment (operating ionizers). (Internal argument set
  automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- downwind_control_expr:

  A quosure (created using
  [`rlang::enquo()`](https://rlang.r-lib.org/reference/enquo.html))
  representing a logical expression used to extract the relevant subset
  of downwind (second stage) observations from `ori_data` that were not
  exposed to treatment (operating ionizers). (Internal argument set
  automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- ori_fitted_models:

  A list containing the models fitted to the `ori_data`, including the
  upwind (first stage) LMM, downwind (second stage) LMM, downwind
  (second stage) treatment-only LMM, downwind (second stage) target-only
  LMM, downwind (second stage) logistic model for rainfall event
  indicator, and downwind (second stage) propensity score model.
  (Internal argument set automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- downwind_lmm_formula:

  A two sided linear formula object to be used in
  [lmer](https://rdrr.io/pkg/lme4/man/lmer.html), describing both the
  fixed-effects and random intercept part of the downwind (second stage)
  LMM. (Internal argument set automatically when using
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

- downwind_propensity_formula:

  A two sided linear formula object to be used in
  [`glm`](https://rdrr.io/r/stats/glm.html) with `family = "binomial"`,
  for fitting a propensity score model to the treatment indicators of
  downwind (second stage) observations. (Internal argument set
  automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- ori_attr_est:

  A numeric vector containing the original attribution estimates (`apo`
  and `apl`) from the original dataset. (Internal argument set
  automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

- ori_sate_est:

  A numeric vector containing the original SATE estimates (`sate.mb`,
  `sate.ipw`, `sate.ipw.l`, `sate.ipw.ma` and `sate.aipw`) from the
  original dataset. (Internal argument set automatically when using
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md))

## Value

A list containing

- hatattr:

  Matrix of bootstrap samples for attribution estimates.

- hatsate:

  Matrix of bootstrap samples for SATE estimates.

- downwind_lmm_param:

  Matrix of bootstrap samples for fixed effect coefficient and random
  effect variance estimates of downwind (second stage) LMM.

- downwind_logistic_param:

  Matrix of bootstrap samples for regression coefficient estimates of
  downwind logistic model fitted to the rainfall event indicators. This
  is `NULL` if `bootstrap_zero = FALSE`.

- downwind_propensity_param:

  Matrix of bootstrap samples for regression coefficient estimates of
  downwind propensity score model fitted to the treatment indicators.
  This is `NULL` if `bootstrap_zero = FALSE`.

- downwind_positive_target_lmm_param:

  Matrix of bootstrap samples for fixed effect coefficient and random
  effect variance estimates of downwind (second stage) treatment-only
  LMM.

- downwind_positive_control_lmm_param:

  Matrix of bootstrap samples for fixed effect coefficient and random
  effect variance estimates of downwind (second stage) control-only LMM.

- downwind_LogRain:

  Matrix of bootstrap samples for the log-transformed rainfall of all
  downwind (second-stage) observations. Observations with zero
  bootstrapped rainfall are represented as `NA`.

## Details

This function implements the bootstrap procedure used in
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).
It is intended for internal use only. Users should not call this
function directly. Instead, bootstrap inference should be performed by
calling
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
with `bootstrap = TRUE`.

**First-Level Bootstrap**  
The first-level is an optional level that is only carried out when
`bootstrap_zero = TRUE`. This level generates bootstrap samples of
binary rainfall event indicator \\L\_{ij}^\*\\ via \\P(L\_{ij}^\* = 1) =
\hat{\gamma}\_{ij} \\ for the subset of observations from `ori_data`
satisfying `downwind`, where \\\hat{\gamma}\_{ij}\\ are the predicted
probabilities from the downwind logistic model fitted to the original
binary rainfall event indicators \\L\_{ij}\\ i.e.,
\\\hat{\gamma}\_{ij}\\ =
`predict(ori_fitted_models$downwind_logistic_fit,type = "response")`. It
is worth noting that if `positive_prob_threshold` is supplied, then
\\\hat{\gamma}\_{ij}\\ that are less than `positive_prob_threshold` are
set to be zeros before being used to generate \\L\_{ij}^\*\\. When
`bootstrap_zero = FALSE`, the first-level bootstrap is not carried out
and thus \\L\_{ij}^\* = L\_{ij}\\.

**Second-Level Bootstrap**  
Recall that the original downwind (second stage) LMM is given as
\$\$y\_{ij} = x\_{ij}^\top \alpha + z\_{ij}^\top \beta + u_i +
e\_{ij},\$\$ the second-level bootstrap generates bootstrap samples of
positive rainfall for the subset of observations from `ori_data`
satisfying `downwind` and \\L\_{ij}^\* = 1\\ from the above first-level
bootstrap via \$\$y\_{ij}^\* = x\_{ij}^\top \hat{\alpha} + z\_{ij}^\top
\hat{\beta} + u_i^\* + e\_{ij}^\* ,\$\$ where \\\hat{\alpha}\\ and
\\\hat{\beta}\\ are estimated fixed effect coefficients from the fitted
downwind (second stage) LMM i.e., `ori_fitted_models$downwind_lmm_fit`,
and \\u_i^\*\\ and \\e\_{ij}^\*\\ are bootstrap samples of the random
intercepts and error terms.

We now explain the procedure for obtaining bootstrap samples of
\\u_i^\*\\ and \\e\_{ij}^\*\\. Let \\\hat{r}\_{ij} = y\_{ij} -
x\_{ij}^\top \hat{\alpha} - z\_{ij}^\top \hat{\beta}\\ be the marginal
residuals for the subset of observations from `ori_data` satisfying
`downwind & ori_positive`, \\\hat{u}\_i = \sum\_{j=1}^{n_i}
\hat{r}\_{ij}/n_i \\ be the day(group)-level average residuals, and
\\\hat{e}\_{ij} = r\_{ij} - \hat{u}\_i \\ be the gauge(unit)-level
residuals, where \\n_i\\ denotes the total number of observations in day
(group) \\i\\ from `ori_data` satisfying `downwind & ori_positive`.
Depending on the chosen `bootstrap_type`, bootstrap samples of
\\u_i^\*\\ and \\e\_{ij}^\*\\ are generated via:

- `REB0`:

  - \\u_i^\* = SRSWR( (\hat{u}\_1, \ldots, \hat{u}\_D ), 1 )\\ for
    \\i=1,\ldots, D^\*\\, where \\D\\ is the number of unique days
    (groups) from `ori_data` satisfying `downwind & ori_positive`,
    \\D^\*\\ is the number of unique days (groups) from `ori_data`
    satisfying `downwind` and \\L\_{ij}^\* = 1\\, and \\SRSWR(a,c)\\
    denote the outcome of \\c\\ independent draws based on simple random
    sampling with replacement from the vector \\a\\.

  - First, sample the donor cluster \\d_i^\* = SRSWR( ( 1,\ldots,D ), 1
    ) \\ for \\i = 1,\ldots, D^\*\\. Then, sample \\e_i^\* =
    (e\_{i1}^\*, \ldots, e\_{in_i^\*}^\*)^\top = SRSWR( (
    \hat{e}\_{d_i^\* 1}, \ldots, \hat{e}\_{d_i^\* n\_{d_i^\*}} ), n_i^\*
    ) \\, where \\n_i^\*\\ denotes the total number of observations in
    day (group) \\i\\ from `ori_data` satisfying `downwind` and
    \\L\_{ij}^\* = 1\\.

- `REB1`:

  Replaces \\\hat{u}\_i\\ and \\\hat{e}\_{ij}\\ in `REB0` with
  \\\hat{u}\_{ij}^{cs} = \hat{\sigma}\_u \hat{u}\_i^{c} \\ D^{-1}
  \sum\_{i' =1}^{D} \hat{u}\_i^2 \\^{-1/2} \\ and \\\hat{e}\_{ij}^{s} =
  \hat{\sigma}\_e \hat{e}\_{ij} \\ N^{-1} \sum\_{i' = 1}^{D} \sum\_{j' =
  1}^{n\_{i'} } \hat{e}\_{i'j'}^2 \\^{-1/2} \\, respectively, where \\
  \hat{u}\_i^c = \hat{u}\_i - D^{-1} \sum\_{i'=1}^{D} \hat{u}\_{i'} \\,
  \\\hat{\sigma}^2_u\\ and \\\hat{\sigma}^2_e\\ are the estimated
  variances of random intercepts and error terms from the fitted
  downwind (second stage) LMM, and \\N = \sum\_{i=1}^{D} n_i\\ is the
  total number of observations from `ori_data` satisfying
  `downwind & ori_positive`.

- `REB2`:

  Same procedure for obtaining \\u_i^\*\\ and \\e\_{ij}^\*\\ as in
  `REB0`, but involves an additional post-processing step on the
  bootstrap estimates of attribution and SATE, as discussed below.

- `PREB0`:

  - \\u_i^\* = SRSWR( (\hat{u}\_1, \ldots, \hat{u}\_D ), 1 )\\ for
    \\i=1,\ldots, D^\*\\.

  - First, sample the donor cluster \\d_i^\* = PPSWR{ (1,\ldots,D),
    (n_1,\cdots, n_D), 1 } \\ for \\i = 1,\ldots, D^\*\\, where
    \\PPSWR(a,b,c)\\ denotes the outcome of \\c\\ independent draws
    based on probability-proportional-to-size sampling with replacement
    from the vector \\a = (a_1,\ldots, a_D)\\ with corresponding sizes
    given by the vector \\b = (b_1,\ldots,b_D)\\, i.e., the probability
    of \\a_i\\ being selected is given as \\b_i / \sum\_{i' = 1}^{D}
    b\_{i'}\\. Then, sample \\e_i^\* = (e\_{i1}^\*, \ldots,
    e\_{in_i^\*}^\*)^\top = SRSWR( ( \hat{e}\_{d_i^\* 1}, \ldots,
    \hat{e}\_{d_i^\* n\_{d_i^\*}} ), n_i^\* ) \\.

- `PREB1`:

  Replaces \\\hat{u}\_i\\ and \\\hat{e}\_{ij}\\ in `PREB0` with
  \\\hat{u}\_{ij}^{sc} = \hat{\sigma}\_u \hat{u}\_i^{c} \\ D^{-1}
  \sum\_{i' =1}^{D} (\hat{u}\_i^c)^2 \\^{-1/2} \\ and
  \\\hat{e}\_{ij}^{s} = \hat{\sigma}\_e \hat{e}\_{ij} \\ N^{-1}
  \sum\_{i' = 1}^{D} \sum\_{j' = 1}^{n\_{i'} } \hat{e}\_{i'j'}^2
  \\^{-1/2} \\, respectively.

- `PREB2`:

  Same procedure for obtaining \\u_i^\*\\ and \\e\_{ij}^\*\\ as in
  `PREB0`, but involves an additional post-processing step on the
  bootstrap estimates of attribution and SATE, as discussed below.

- `MREB1`:

  Replaces \\\hat{u}\_i\\ and \\\hat{e}\_{ij}\\ in `REB0` with
  \\\hat{u}\_{ij}^{sc}\\ defined under `PREB1` and \\\tilde{e}\_{ij}^{s}
  = \hat{\sigma}\_e \hat{e}\_{ij} \\ \sum\_{i' = 1}^{D} \sum\_{j' =
  1}^{n\_{i'} } D^{-1} n\_{i'}^{-1} \hat{e}\_{i'j'}^2 \\^{-1/2} \\,
  respectively.

The random effect block (REB0, REB1, REB2) bootstraps proposed in
Chambers and Chandra (2013) were originally designed to handle balanced
clustered data, while the proportional REB (PREB0, PREB1, PREB2)
bootstraps and the MREB1 bootstrap proposed by Tho et al. (2025) are
generalizations of the REB bootstraps to accommodate highly unbalanced
clustered data. Therefore, it is recommended to use either
`bootstrap_type = "PREB1"` or `bootstrap_type = "MREB1"`, especially
when \\n_i\\'s are highly unbalanced. Users are refered to Tho et al.
(2025) for more discussion on the comparison among these bootstrap
methods.

**Adjustment of Bootstrapped Rainfall**  
After obtaining the bootstrapped \\y\_{ij}^\*\\ that are assumed to be
log-rainfall, it is possible to perform some adjustment to the bootstrap
samples to ensure that they are similar to the observed data.

Let \\Rain\_{ij}^\* = \exp(y\_{ij}^\*)\\ be the bootstrapped raw
rainfall. When `discretize_rain = TRUE`, bootstrap samples that satisfy
\\ Rain\_{ij}^\* \in (0, 0.3\] \\ are replaced by 0.2, \\ Rain\_{ij}^\*
\in (0.3, 0.5\] \\ are replaced by 0.4, \\ Rain\_{ij}^\* \in (0.5, 0.7\]
\\ are replaced by 0.6, and \\ Rain\_{ij}^\* \in (0.7, 0.9\] \\ are
replaced by 0.8.

When `winsorize_individual_rain = TRUE`, bootstrap samples that satisfy
\\Rain\_{ij}^\* \> \\ `individual_rain_interval[2]` are replaced by
random numbers drawn from a uniform distribution over the interval
\\\[\\`individual_rain_interval[1]`,
`individual_rain_interval[2]`\\\]\\.

When `winsorize_total_rain = TRUE`, if \\ \sum\_{(i,j)} Rain\_{ij}^\*
\notin \[\\ `total_rain_interval[1]`, `total_rain_interval[2]` \\\]\\
where the summation is over the subset of observations in `ori_data`
satisfying `downwind` and \\L\_{ij}^\* = 1\\, then all bootstrap samples
of \\Rain\_{ij}^\*\\ are rescaled by a common factor of
`runif(n = 1, min = total_rain_interval[1], max = total_rain_interval[2])`
\\ / \sum\_{(i,j)} Rain\_{ij}^\* \\.

The final adjusted \\Rain\_{ij}^\*\\ are converted back to the log-scale
based on the formula \\y\_{ij}^\* = \log(Rain\_{ij}^\*)\\. Therefore,
these adjustments should only be used when modelling log-transformed
rainfall in the two-stage LMM approach, but not raw rainfall. If an
offset term is included on the LHS of `downwind_lmm_formula` e.g.,
`downwind_lmm_formula = LogRain - Offset`, the above adjustments could
still be used as the function would use the relationship \\Rain\_{ij}^\*
= \exp( y\_{ij}^\* + Offset\_{ij} )\\ and \\y\_{ij}^\* =
\log(Rain\_{ij}^\*) + Offset\_{ij} \\.

**Bootstrap Distributions of Attribution, SATE and other parameters**  
The same procedure described in
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
is then repeated on the bootstrapped dataset, where the bootstrap
samples of \\y\_{ij}^\*\\ generated from the above two-level bootstrap
procedure are used to replace the original observed \\y\_{ij}\\. This
includes the fitting of the downwind (second stage) LMM, downwind
(second-stage) target-only LMM, and downwind (second-stage) control-only
LMM to the subset of observations from `ori_data` satisfying
`downwind &` \\ L\_{ij}^\* = 1\\. When `bootstrap_zero = TRUE`, two
additional models are also fitted; namely, the downwind logistic model
fitted to the subset of observations from `ori_data` satisfying
`downwind` with the response being the bootstrapped rainfall event
indicator \\L\_{ij}^\*\\, and the downwind propensity score model fitted
to the subset of observations from `ori_data` satisfying `downwind &`
\\L\_{ij}^\* = 1\\ with the response being the indicator \\I\_{ij}\\ for
exposure to ionizer (treatment). Finally, two attribution estimates and
five SATE estimates are computed based on the estimation results of
these models fitted to the bootstrapped dataset, where \\Rain\_{ij}\\
and \\y\_{ij}\\ are replaced by their bootstrap counterparts
\\Rain\_{ij}^\*\\ and \\y\_{ij}^\*\\, respectively.

By repeatedly generating bootstrap samples, fitting models and computing
attribution and SATE estimates for `B_bootstrap` number of times, this
function returns the bootstrap distributions of

- Two attribution estimates

- Five SATE estimates

- Fixed effect coefficients and random effect variance estimates of
  downwind LMM, downwind treatment-only LMM, and downwind control-only
  LMM.

- Log-transformed rainfall \\\log(Rain\_{ij}^\*)\\

- Regression coefficient estimates of downwind logistic model fitted to
  the bootstrapped rainfall event indicators, only when
  `bootstrap_zero = TRUE`

- Regression coefficient estimates of downwind propensity score model,
  only when `bootstrap_zero = TRUE`

**Parallel Bootstrap Execution**  
This function also supports parallel execution via the argument
`bootstrap_parallel`. When `bootstrap_parallel = TRUE`, the
`B_bootstrap` repeated iterations are executed concurrently across
`bootstrap_parallel_num_worker` workers. Each iteration, including
rainfall resampling, model fitting, and parameter estimation, is
performed independently on different workers, which can reduce
computation time compared to sequential execution. The parallel
execution is implemented using the
[`foreach`](https://rdrr.io/pkg/foreach/man/foreach.html) package with
the
[`doParallel`](https://rdrr.io/pkg/doParallel/man/doParallel-package.html)
backend and
[`registerDoRNG`](https://rdrr.io/pkg/doRNG/man/registerDoRNG.html) for
reproducibility.

**Post-processing of `REB2` and `PREB2`**  
The `REB2` and `PREB2` perform post-processing on all of the above
bootstrap distributions (except for log-transformed rainfall) generated
by `REB0` and `PREB0`, respectively.

Specifically:

- The bootstrap distributions of attribution, SATE, fixed-effect
  coefficients from LMMs, and regression coefficients from the downwind
  logistic model and downwind propensity score model are mean-centered
  at their original estimates using: \$\$\hat{\theta}^\*\_b \leftarrow
  \hat{\theta}^\*\_b + \hat{\theta} -
  \frac{1}{B}\sum\_{b=1}^{B}\hat{\theta}^\*\_b\$\$ where
  \\\hat{\theta}^\*\_b\\ is the estimate from the \\b\\-th bootstrap
  sample and \\\hat{\theta}\\ is the corresponding estimate from the
  original dataset.

- For LMM random effect variance components, the bootstrap distributions
  of random intercept variance (\\\sigma^2_u\\) and residual variance
  (\\\sigma^2_e\\) are first adjusted to be empirically uncorrelated
  (see Section 2.3.3 of Chambers & Chandra (2013), and Section 2.3 of
  Tho et al. (2025)).

- After adjustment, ratio corrections are applied:
  \$\$\hat{\sigma}^{2\*}\_{u,b} \leftarrow \hat{\sigma}^{2\*}\_{u,b}
  \times
  \frac{\hat{\sigma}^2_u}{\frac{1}{B}\sum\_{b=1}^{B}\hat{\sigma}^{2\*}\_{u,b}}\$\$
  \$\$\hat{\sigma}^{2\*}\_{e,b} \leftarrow \hat{\sigma}^{2\*}\_{e,b}
  \times
  \frac{\hat{\sigma}^2_e}{\frac{1}{B}\sum\_{b=1}^{B}\hat{\sigma}^{2\*}\_{e,b}}\$\$
  where \\\hat{\sigma}^{2\*}\_{u,b}\\ and \\\hat{\sigma}^{2\*}\_{e,b}\\
  are the adjusted bootstrap estimates, and \\\hat{\sigma}^2_u\\,
  \\\hat{\sigma}^2_e\\ are the original estimates.

**Bootstrap P-values, Bootstrap Percentile Confidence Intervals,
Bootstrap Plots**  
All bootstrap distributions (except for log-transformed rainfall)
produced by this function are further used in
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
to:

- Compute bootstrap p-values as the proportion of bootstrapped estimates
  that are less than zero, i.e., \$\$ \frac{1}{B} \sum\_{b=1}^{B} 1\_{
  \\ \hat{\theta}\_b^\* \< 0 \\} \$\$ where \\1\_{\\\cdot\\}\\ is the
  indicator function and \\\hat{\theta}\_b\\ denotes the estimate from
  the \\b\\-th bootstrap sample.

- Compute bootstrap percentile confidence interval with confidence level
  \\(1-\alpha) \times 100\\\\ as \$\$ \[ \hat{\theta}^\*\_{\alpha / 2},
  \hat{\theta}^\*\_{1 - \alpha/2} \], \$\$ where
  \\\hat{\theta}^\*\_{p}\\ is the empirical \\p\\-th quantile of the
  bootstrap distribution \\\\ \hat{\theta}^\*\_{1}, \ldots,
  \hat{\theta}^\*\_{B} \\\\ of the parameter estimate \\\hat{\theta}\\

Additionally, the bootstrap distributions of attribution and SATE are
also plotted using kernel density estimates, with:

- A solid vertical line for the original estimate \\\hat{\theta}\\;

- A dotted vertical line at zero for reference.

## References

- Chambers, R. and Chandra, H. (2013). A Random Effect Block Bootstrap
  for Clustered Data. *Journal of Computational and Graphical
  Statistics*, 22, 452–470.

- Tho, Z. Y., Chambers, R., and Welsh, A. H. (2025) A Proportional
  Random Effect Block Bootstrap for General Clustered Data.
  [arXiv:2510.07770](https://arxiv.org/abs/2510.07770).

## See also

[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
for the main function,
[`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md)
for specifying bootstrap options
