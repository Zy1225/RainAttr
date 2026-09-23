# Attribution and Sample Average Treatment Effect for Rainfall Enhancement Trial Data

Perform estimation and inference of attribution and sample average
treatment effect for rainfall enhancement trial data, based on the
two-stage linear mixed model (LMM) approach employed in Chambers et al.
(2022a).

## Usage

``` r
rain_attr(
  data,
  upwind_lmm_formula,
  instr_pred_name,
  instr_pred_type,
  downwind_lmm_formula,
  downwind_logistic_formula = NULL,
  downwind_propensity_formula,
  rain_col_name,
  upwind_subset,
  downwind_subset,
  downwind_target_subset,
  downwind_control_subset,
  positive_subset,
  attr_type = "ThoEtAl",
  x_downwind_name,
  target_only = FALSE,
  bootstrap = FALSE,
  bootstrap_option = bootstrap_opt(),
  permutation = FALSE,
  permutation_option = permutation_opt()
)
```

## Arguments

- data:

  A data frame containing the variables named in `upwind_lmm_formula`,
  `downwind_lmm_formula`, `downwind_logistic_formula` (if specified),
  and `downwind_propensity_formula`. It should also contain variables
  named in `rain_col_name`, `upwind_subset`, `downwind_subset`,
  `downwind_target_subset`, and `downwind_control_subset`.

- upwind_lmm_formula:

  A two sided linear formula object to be used in
  [lmer](https://rdrr.io/pkg/lme4/man/lmer.html), describing both the
  fixed-effects and random intercept part of the upwind (first stage)
  LMM.

- instr_pred_name:

  A character string to store the variable name of the fitted values
  generated from the upwind (first stage) LMM.

- instr_pred_type:

  Type of fitted values generated from the upwind (first stage) LMM. If
  "Unconditional" the fitted values equal to only the estimated fixed
  effects. If "Conditional" the fitted values equal to the sum of
  estimated fixed effects and EBLUPs of random intercepts.

- downwind_lmm_formula:

  A two sided linear formula object to be used in
  [lmer](https://rdrr.io/pkg/lme4/man/lmer.html), describing both the
  fixed-effects and random intercept part of the downwind (second stage)
  LMM. This formula should contain the variable name specified in
  `instr_pred_name`.

- downwind_logistic_formula:

  An optional two sided linear formula object to be used in
  [`glm`](https://rdrr.io/r/stats/glm.html) with `family = "binomial"`,
  for fitting a logistic model to the indicators of rainfall event. This
  only needs to be specified when `bootstrap = TRUE` and
  `bootstrap_option$bootstrap_zero = TRUE`.

- downwind_propensity_formula:

  A two sided linear formula object to be used in
  [`glm`](https://rdrr.io/r/stats/glm.html) with `family = "binomial"`,
  for fitting a propensity score model to the treatment indicators of
  downwind (second stage) observations.

- rain_col_name:

  A character string that refers to the column name of the raw scale
  rainfall in `data`.

- upwind_subset:

  A logical expression used to extract the relevant subset of
  observations from `data` to be used in the upwind (first stage) LMM
  fitting. For example, `Gauge.Day.Type == "Upwind"`.

- downwind_subset:

  A logical expression used to extract the relevant subset of
  observations from `data` to be used in the downwind (second stage) LMM
  fitting. For example, `Gauge.Day.Type %in% c("Target","Control")`.

- downwind_target_subset:

  A logical expression used to extract the relevant subset of downwind
  (second stage) observations from `data` that were exposed to treatment
  (operating ionizers). For example, `Gauge.Day.Type == "Target"`.

- downwind_control_subset:

  A logical expression used to extract the relevant subset of downwind
  (second stage) observations from `data` that were not exposed to
  treatment (operating ionizers). For example,
  `Gauge.Day.Type == "Control"`.

- positive_subset:

  A logical expression used to extract the relevant subset of
  observations from `data` with positive rainfall - these are the
  observations that are used in the fitting of upwind (first stage) LMM,
  downwind (second stage) LMM, downwind (second stage) treatment-only
  LMM, downwind (second stage) control-only LMM, and the downwind
  (second stage) propensity score model.

- attr_type:

  An optional character string specifying the type of attribution
  estimates. Must be one of `"ChambersEtAl"`,
  `"ChambersEtAl_No_Winsorize"`, `"ThoEtAl"` (default), or `"No"`. See
  "Details" for more information.

- x_downwind_name:

  A character vector containing variable names from the right hand side
  of `downwind_lmm_formula`, for those variables that are not related to
  ionizers (treatment). The intercept is always included and does not
  need to be specified.

- target_only:

  An optional logical. If `TRUE` the attribution estimates are computed
  based on only treated observations. If `FALSE` the attribution
  estimates are computed based on both treated and control observations.

- bootstrap:

  An optional logical. If `TRUE` bootstrap is carried out to perform
  inference on the attribution and sample average treatment effect. If
  `FALSE` (default) no bootstrap is carried out.

- bootstrap_option:

  An optional list containing all bootstrap settings, used only when
  `bootstrap = TRUE`. See
  [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md)
  for the default list elements and their usage.

- permutation:

  An optional logical, If `TRUE` randomized permutation is carried out
  on the ionizer operation (treatment) schedule to perform inference on
  the attribution and sample average treatment effect. If `FALSE`
  (default) no randomized permutation is carried out.

- permutation_option:

  An optional list containing all permutation settings, used only when
  `permutation = TRUE`. See
  [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md)
  for the default list elements and their usage.

## Value

An object of class
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md),
which is a list containing

- all_fitted_models:

  A list of model objects from
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) for the first stage
  (upwind), second stage (downwind) LMM, second stage (downwind)
  treatment-only LMM, and second stage (downwind) control-only LMM,
  along with model objects from
  [`glm`](https://rdrr.io/r/stats/glm.html) for the logistic model of
  rainfall event indicator (`NULL` if `downwind_logistic_formula` is not
  specified) and the propensity score model for the treatment indicator
  of second stage (downwind) observations.

- hatattr:

  A vector containing the attribution estimates.

- hatsate:

  A vector containing the sample average treatment effect estimates.

- bootstrap_result:

  A list of matrices with the following elements:

  - hatattr: Matrix of bootstrap samples for attribution estimates.

  - hatsate: Matrix of bootstrap samples for SATE estimates.

  - downwind_lmm_param: Matrix of bootstrap samples for fixed effect
    coefficient and random effect variance estimates of downwind (second
    stage) LMM.

  - downwind_logistic_param: Matrix of bootstrap samples for regression
    coefficient estimates of downwind logistic model fitted to the
    rainfall event indicators. This is `NULL` if
    `downwind_logistic_formula` is not specified.

  - downwind_propensity_param: Matrix of bootstrap samples for
    regression coefficient estimates of downwind propensity score model
    fitted to the treatment indicators.

  - downwind_positive_target_lmm_param: Matrix of bootstrap samples for
    fixed effect coefficient and random effect variance estimates of
    downwind (second stage) treatment-only LMM.

  - downwind_positive_control_lmm_param: Matrix of bootstrap samples for
    fixed effect coefficient and random effect variance estimates of
    downwind (second stage) control-only LMM.

  - downwind_LogRain: Matrix of bootstrap samples for the
    log-transformed rainfall of all downwind (second-stage)
    observations. Observations with zero bootstrapped rainfall are
    represented as `NA`.

  This is `NULL` when `bootstrap = FALSE`.

- bootstrap_CI_result:

  A list of matrices with same element names as in `bootstrap_result`
  (excluding `downwind_LogRain`), containing the corresponding bootstrap
  percentile confidence intervals. This is `NULL` when
  `bootstrap = FALSE`.

- bootstrap_p_value_result:

  A list of numeric vectors with same element names as in
  `bootstrap_result` (excluding `downwind_LogRain`), containing the
  corresponding proportion of bootstrap samples that are less than zero.
  This is `NULL` when `bootstrap = FALSE`.

- bootstrap_plot_result:

  A list with two elements:

  - hatattr: A list of `ggplot` objects, each showing the bootstrap
    distribution of attribution estimates. Each plot includes a dotted
    vertical line at zero and a solid vertical line at the original
    estimate based on the observed data.

  - hatsate: A list of `ggplot` objects, each showing the bootstrap
    distribution of SATE estimates. Each plot includes a dotted vertical
    line at zero and a solid vertical line at the original estimate
    based on the observed data. This is `NULL` when `bootstrap = FALSE`.

- permutation_result:

  A list of matrices with two elements:

  - hatattr: Matrix of permutation samples for attribution estimates.

  - hatsate: Matrix of permutation samples for SATE estimates.

  This is `NULL` when `permutation = FALSE`.

- permutation_p_value_result:

  A list of numeric vectors with same element names as in
  `permutation_result`, containing the corresponding proportion of
  permutation samples that are greater than or equal to the original
  estimate based on the observed data. This is `NULL` when
  `permutation = FALSE`.

- permutation_plot_result:

  A list with two elements:

  - hatattr: A list of `ggplot` objects, each showing the permutation
    distribution of attribution estimates. Each plot includes a solid
    vertical line at the original estimate based on the observed data.

  - hatsate: A list of `ggplot` objects, each showing the permutation
    distribution of SATE estimates. Each plot includes a solid vertical
    line at the original estimate based on the observed data.

  This is `NULL` when `permutation = FALSE`.

- args:

  A list of the original function arguments.

- data:

  A data frame containing the original supplied `data`, with an
  additional column containing the fitted values generated from the
  upwind (first stage) LMM.

## Details

This function implements a two-stage modelling procedure via the
following steps:

1.  Fit an upwind (first stage) LMM using
    `lme4::lmer(upwind_lmm_formula)` to the subset of observations from
    `data` satisfying `upwind_subset & positive_subset`. This fitted LMM
    is used to obtain fitted values (named as `instr_pred_name`) to be
    used in the downwind (second stage) LMM.

2.  Fit a downwind (second stage) LMM using
    `lme4::lmer(downwind_lmm_formula)` to the subset of observations
    from `data` satisfying `downwind_subset & positive_subset`.
    \$\$y\_{ij} = x\_{ij}^\top \alpha + z\_{ij}^\top \beta + u_i +
    e\_{ij} \$\$ where \\i\\ indexes day (group) and \\j\\ indexes gauge
    (unit within group), \\x\_{ij}\\ is a vector of covariates (with
    names supplied in `x_downwind_name`, including intercept) that are
    not related to the ionizers (treatment), \\z\_{ij}\\ is a vector of
    ionizer (treatment) related covariates, \\u_i\\ are random
    intercepts, and \\e\_{ij}\\ are error terms.

The fitted values obtained from the upwind (first stage) LMM can either
be:

- Included as a covariate on the right-hand side of
  `downwind_lmm_formula`, e.g., `instr_pred_name = "natural_pred"` and
  `downwind_lmm_formula = LogRain ~ natural_pred + ...`, or

- Included as an offset term by subtracting it from the response on the
  left-hand side, e.g., `instr_pred_name = "natural_pred"` and
  `downwind_lmm_formula = LogRain - natural_pred ~ ...`.

**Attribution**  
Two attribution estimates, namely `apo` and `apl` are computed based on
the estimated fixed effect coefficients \\\hat{\alpha}\\,
\\\hat{\beta}\\ and EBLUPs \\\hat{u}\_i\\ from the fitted downwind
(second stage) LMM. `apo` represents the total increase or decrease in
downwind rainfall attributed to the ionizer (treatment) as a proportion
of the total amount of observed downwind rainfall, while `apl`
represents the total increase or decrease in downwind rainfall
attributed to the ionizer (treatment) as a proportion of the total
expected amount of downwind rainfall without the effect of ionizer
(treatment). This function allows for four different ways of estimating
`apo` and `apl` as specified by the argument `attr_type`:

- `ChambersEtAl`:

  Attribution is estimated based on the approach of Chambers et al.
  (2022a), to adjust for back-transformation bias due to the modelling
  of log-transformed rainfall: \$\$ \code{apo} = \sum\_{(i,j)}
  Rain\_{ij} \[ 1 - \max\\\lambda^{-1} \exp(-z\_{ij}^\top \hat{\beta}),
  0.5\\ \] / \sum\_{(i,j)}Rain\_{ij}, \$\$ \$\$ \code{apl} =
  \sum\_{(i,j)} Rain\_{ij} \[ 1 - \max\\\lambda^{-1} \exp(-z\_{ij}^\top
  \hat{\beta}), 0.5\\ \] / \sum\_{(i,j)}Rain\_{ij} \max\\\lambda^{-1}
  \exp(-z\_{ij}^\top \hat{\beta}), 0.5\\, \$\$ where the summation is
  either across all observations satisfying
  `downwind_subset & positive_subset` (when `target_only = FALSE`), or
  across all observations satisfying
  `downwind_target_subset & positive_subset` (when
  `target_only = TRUE`), \\Rain\_{ij}\\ is the observed raw rainfall
  (contained in the column specified by `rain_col_name`), \$\$ \lambda =
  1 + \frac{\sqrt{ (1+m)^2 + 4(\mu - 1)m } - (1+m)}{2m}, m =
  \frac{\hat{V}( x\_{ij}^\top \hat{\alpha} + \hat{u}\_i )
  }{\hat{V}(z\_{ij}^\top \hat{\beta})}, \$\$ with \\\hat{V}(\cdot)\\
  denoting the empirical variance either across all observations
  satisfying `downwind_subset & positive_subset` (when
  `target_only = FALSE`), or across all observations satisfying
  `downwind_target_subset & positive_subset` (when
  `target_only = TRUE`), and \$\$ \mu = \frac{1}{N} \sum\_{(i,j)}
  \frac{Rain\_{ij}}{\exp( x\_{ij}^\top \hat{\alpha} + z\_{ij}^\top
  \hat{\beta} + \hat{u}\_i )}, \$\$ and \\N\\ is either the total number
  of observations satisfying `downwind_subset & positive_subset` (when
  `target_only = FALSE`), or the total number of observations satisfying
  `downwind_target_subset & positive_subset` (when
  `target_only = TRUE`). When an offset term is included on the LHS of
  `downwind_lmm_formula`, the expressions of \\m\\ and \\\mu\\ become
  \$\$ m = \frac{\hat{V}( offset\_{ij} + x\_{ij}^\top \hat{\alpha} +
  \hat{u}\_i ) }{\hat{V}(z\_{ij}^\top \hat{\beta})}, \mu = \frac{1}{N}
  \sum\_{(i,j)} \frac{Rain\_{ij}}{\exp( offset\_{ij} + x\_{ij}^\top
  \hat{\alpha} + z\_{ij}^\top \hat{\beta} + \hat{u}\_i )}. \$\$

- `ChambersEtAl_No_Winsorize`:

  Similar to `ChambersEtAl`, but without the winsorizing step, i.e.,
  \\\max\\\lambda^{-1} \exp(-z\_{ij}^\top \hat{\beta}), 0.5\\\\ is
  replaced by \\\lambda^{-1} \exp(-z\_{ij}^\top \hat{\beta})\\: \$\$
  \code{apo} = \sum\_{(i,j)} Rain\_{ij} \\ 1 - \lambda^{-1}
  \exp(-z\_{ij}^\top \hat{\beta}) \\ / \sum\_{(i,j)}Rain\_{ij}, \quad
  \$\$ \$\$ \code{apl} = \sum\_{(i,j)} Rain\_{ij} \\ 1 - \lambda^{-1}
  \exp(-z\_{ij}^\top \hat{\beta}) \\ / \sum\_{(i,j)}Rain\_{ij}
  \lambda^{-1} \exp(-z\_{ij}^\top \hat{\beta}), \$\$ where \\\lambda\\
  is the same as in `ChambersEtAl`.

- `ThoEtAl`:

  Attribution is estimated based on an alternative adjustment proposed
  by Tho et al. (2026), using the estimated covariance matrix
  \\\hat{\Sigma}\\ of \\\hat{\beta}\\. \$\$ \code{apo} = \sum\_{(i,j)}
  Rain\_{ij} \\ 1 - \exp(-z\_{ij}^\top \hat{\beta} - 0.5 z\_{ij}^\top
  \hat{\Sigma} z\_{ij} ) \\/ \sum\_{(i,j)}Rain\_{ij}, \quad \$\$ \$\$
  \code{apl} = \sum\_{(i,j)} Rain\_{ij} \\ 1 - \exp(-z\_{ij}^\top
  \hat{\beta} - 0.5 z\_{ij}^\top \hat{\Sigma} z\_{ij} ) \\ /
  \sum\_{(i,j)}Rain\_{ij} \exp(-z\_{ij}^\top \hat{\beta} - 0.5
  z\_{ij}^\top \hat{\Sigma} z\_{ij} ). \$\$

- `No`:

  Attribution is estimated based on no adjustment. \$\$ \code{apo} =
  \sum\_{(i,j)} Rain\_{ij} \\ 1 - \exp(-z\_{ij}^\top \hat{\beta} ) \\/
  \sum\_{(i,j)}Rain\_{ij}, \quad \$\$ \$\$\code{apl} = \sum\_{(i,j)}
  Rain\_{ij} \\ 1 - \exp(-z\_{ij}^\top \hat{\beta} ) \\ /
  \sum\_{(i,j)}Rain\_{ij} \exp(-z\_{ij}^\top \hat{\beta} ). \$\$

**SATE**  
The computation of SATE estimates involves fitting a downwind (second
stage) propensity score model using
`glm(downwind_propensity_formula, family = "binomial")` to the subset of
observations from `data` satisfying `downwind_subset & positive_subset`,
with the response being an indicator \\I\_{ij}\\ for whether each
observation is exposed to the ionizer (treatment), i.e., \\I\_{ij} = 1\\
if it satisfies `downwind_target_subset & positive_subset`, and
\\I\_{ij} = 0\\ if it satisfies
`downwind_control_subset & positive_subset`. The estimated propensity
scores (i.e., fitted values) from this fitted propensity score model,
denoted as \\\hat{\pi}\_{ij}\\, are then used to compute the inverse
propensity weights (IPW) \\\hat{w}\_{ij,1} = \hat{\pi}\_{ij}^{-1} /
\sum\_{(k,l)} (\hat{\pi}\_{kl}^{-1} I\_{kl}) \\ and \\ \hat{w}\_{ij,0} =
(1 - \hat{\pi}\_{ij})^{-1} / \sum\_{(k,l)} \\ (1- \hat{\pi}\_{kl})^{-1}
(1- I\_{kl}) \\ \\, where the summation is over all observations from
`data` satisfying `downwind_subset & positive_subset`. These IPW weights
are then used, together with the estimation results of the downwind
(second stage) LMM, to obtain the following five types of SATE estimates
discussed in Chambers et al. (2022b):

- `sate.mb` \\ = \sum\_{(i,j)} I\_{ij} (z\_{ij}^\top \hat{\beta}) /
  \sum\_{(i,j)} I\_{ij} \\.

- `sate.ipw` \\ = \\ \sum\_{(i,j)} \hat{w}\_{ij,1} I\_{ij} y\_{ij} \\ -
  \\ \sum\_{(i,j)} \hat{w}\_{ij,0} (1-I\_{ij}) y\_{ij} \\ \\, where
  \\y\_{ij}\\ denote the response variable of the downwind (second
  stage) LMM that might contain offset term.

- `sate.ipw.l` \\ = \sum\_{(i,j)} \hat{w}\_{ij,1} I\_{ij} z\_{ij}^\top
  \hat{\beta} \\.

- `sate.ipw.ma` \\ = (N^{-1} \sum\_{(i,j)} \hat{m}\_{ij,1} ) - (N^{-1}
  \sum\_{(i,j)} \hat{m}\_{ij,0} ) + \\ \sum\_{(i,j)} \hat{w}\_{ij,1}
  I\_{ij} (y\_{ij} - \hat{m}\_{ij,1}) \\ - \\ \sum\_{(i,j)}
  \hat{w}\_{ij,0} (1- I\_{ij}) (y\_{ij} - \hat{m}\_{ij,0}) \\ \\, where
  \\N\\ is the total number of observations satisfying
  `downwind_subset & positive_subset`. \\\hat{m}\_{ij,1}\\ are fitted
  values (using fixed effect only) obtained from the downwind (second
  stage) treatment-only LMM, based on a modified version of
  `downwind_lmm_formula` whose RHS only contains non-ionizer
  (non-treatment) related covariates \\x\_{ij}\\ and excludes ionizer
  (treatment) related covariates \\z\_{ij}\\ fitted to the subset of
  observations from `data` satisfying
  `downwind_target_subset & positive_subset`. Similarly,
  \\\hat{m}\_{ij,0}\\ are fitted values (using fixed effect only)
  obtained from the downwind (second stage) control-only LMM, based on
  the same modified version of `downwind_lmm_formula` fitted to another
  subset of observations from `data` satisfying
  `downwind_control_subset & positive_subset`.

- `sate.aipw` \\ = \[\sum\_{(i,j)} \hat{w}\_{ij,1} \\ I\_{ij} y\_{ij} -
  (I\_{ij} - \hat{\pi}\_{ij}) \hat{m}\_{ij,1} \\ \] - \[\sum\_{(i,j)}
  \hat{w}\_{ij,0} \\ (1 - I\_{ij}) y\_{ij} - (I\_{ij} - \hat{\pi}\_{ij})
  \hat{m}\_{ij,0} \\ \] \\.

**Bootstrap and Permutation Inference**  
This function can also be used to perform bootstrap inference on the
attribution and SATE, by setting `bootstrap = TRUE` and supplying the
relevant bootstrap options using
[`bootstrap_opt()`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md).
For full details of the bootstrap procedure, please see
[`bootstrap_downwind`](https://zy1225.github.io/RainAttr/reference/bootstrap_downwind.md).
Briefly, the bootstrap is carried out in two levels by conditioning on
the upwind (first stage) LMM and its fitted values:

- First level is an optional level that is only carried out when
  `bootstrap_opt(bootstrap_zero = TRUE)`. This level considers
  generating bootstrap samples of rainfall event indicator for the
  subset of observations satisfying `downwind_subset` using the
  predicted probabilities from the downwind logistic model. This model
  is fitted using `glm(downwind_logistic_formula, family = "binomial")`
  to the subset of observations from `data` satisfying
  `downwind_subset`, with the response being an indicator for whether
  the observed rainfall is greater than zero, i.e., the indicator is
  defined by the logical expression `positive_subset`.

- Second level generates bootstrap samples of positive rainfall for the
  subset of observations not only satisfying `downwind_subset` but also
  with the first-level bootstrapped rainfall event indicator being equal
  to one. When `bootstrap_opt(bootstrap_zero = FALSE)`, then this level
  generates bootstrap samples of positive rainfall for the subset of
  observations satisfying `downwind_subset & positive_subset`. This is
  done using one of the semiparametric bootstrap methods of Chambers &
  Chandra (2013) and Tho et al. (2025), which involves the use of
  marginal residuals from the fitted downwind (second stage) LMM.

The above attribution and SATE estimates are then computed based on each
bootstrap sample of the positive rainfall, forming their respective
bootstrap distributions. This function also provides bootstrap
distributions of parameters associated with the downwind LMM
(`downwind_lmm_formula`), downwind logistic model
(`downwind_logistic_formula`), downwind propensity score model
(`downwind_propensity_formula`), downwind treatment-only LMM, and
downwind control-only LMM. These bootstrap distributions are then used
to compute bootstrap p-values (proportion of bootstrapped estimates that
are less than or equal to zero), form bootstrap percentile confidence
intervals (with confidence level specified in
`bootstrap_option$CI_level`), and generate their respective plots. It is
worth noting that the entire bootstrap procedure (including rainfall
resampling, model fitting, and parameter estimation) can be run in
parallel by setting `bootstrap_option$bootstrap_parallel = TRUE`, using
`bootstrap_option$bootstrap_parallel_num_worker` workers.

Finally, this function enables permutation-based inference on the
attribution and SATE, by setting `permutation = TRUE` and supplying the
relevant permutation options using
[`permutation_opt()`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md).
For full details of the permutation-based procedure, please see
[`permutation_ionizer`](https://zy1225.github.io/RainAttr/reference/permutation_ionizer.md).
In short, the permutation-based procedure involves randomly permuting
the operating schedules of the ionizers (treatment) and re-estimating
the attribution and SATE based on the permuted data, from which
permutation distributions of attribution and SATE estimates are formed.
These permutation distributions are used to compute permutation p-values
(proportion of permuted estimates that are greater than or equal to the
observed estimates) and generate their respective plots.

## References

- Chambers, R., Beare, S., Peak, S. and Al-Kalbani, M. (2022a) Nudging a
  Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall
  Enhancement Trial in Oman. *International Statistical Review*, 90:
  346–373.

- Chambers, R. and Chandra, H. (2013). A Random Effect Block Bootstrap
  for Clustered Data. *Journal of Computational and Graphical
  Statistics*, 22, 452–470.

- Chambers, R., Ranjbar, S., Salvati, N., and Pacini, B. (2022b)
  Weighting, Informativeness and Causal Inference, with an Application
  to Rainfall Enhancement. *Journal of the Royal Statistical Society
  Series A: Statistics in Society*, 185: 1584–1612

- Tho, Z. Y., Chambers, R., and Welsh, A. H. (2025) A Proportional
  Random Effect Block Bootstrap for General Clustered Data.
  [arXiv:2510.07770](https://arxiv.org/abs/2510.07770).

- Tho, Z. Y., Chambers, R., and Welsh, A. H. (2026) Bias-Adjusted
  Attribution Estimation for Rainfall Enhancement Trials.

## Examples

``` r
result <- rain_attr(
  data = oman,
  upwind_lmm_formula =
    LogRain ~ Gauge.Elevation + Steering.Wind.Speed + Total.Totals +
    PC2.Dry.Temperature + PC1.Relative.Humidity +
    PC1.Ground.Level.Pressure + (1 | TrialDay),
  instr_pred_name = "natural_pred",
  instr_pred_type = "Conditional",
  downwind_lmm_formula =
    LogRain - natural_pred ~ Gauge.Elevation +
    Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 +
    Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 +
    Target.H.09 + Target.H.10 +
    Gauge.Elevation:Target.H.01 +
    Gauge.Elevation:Target.H.02 +
    (1 | TrialDay),
  downwind_logistic_formula = NULL,
  downwind_propensity_formula =
    (Gauge.Day.Type == "Target") ~ Total.Totals +
    PC1.Dry.Temperature + PC1.Relative.Humidity +
    PC1.Ground.Level.Pressure,
  rain_col_name = "Rain.Gauge.Measurement",
  upwind_subset = Gauge.Day.Type == "Upwind",
  downwind_subset = Gauge.Day.Type %in% c("Target", "Control"),
  downwind_target_subset = Gauge.Day.Type == "Target",
  downwind_control_subset = Gauge.Day.Type == "Control",
  positive_subset = Rain.Gauge.Measurement > 0,
  attr_type = "ThoEtAl",
  x_downwind_name = "Gauge.Elevation",
  bootstrap = FALSE,
  permutation = FALSE
)

# Print and summarize the fitted analysis
print(result)
#> Two Stage LMM Rainfall Enhancement Analysis Result
#> ======================================================================
#> 
#> Point Estimates:
#> Attribution (%) Assuming Log-Rainfall being Modelled:
#>   apo   apl 
#> 6.27% 6.69% 
#> 
#> SATE Estimates (hatsate):
#>     sate.mb    sate.ipw  sate.ipw.l sate.ipw.ma   sate.aipw 
#>     0.11531     0.07513     0.11324     0.08135     0.07706 
#> ======================================================================
#> 
#> Inference:
#> Bootstrap inference has NOT been carried out.
#> 
#> Permutation inference has NOT been carried out.
#> 
#> ======================================================================
#> 
#> Upwind (First Stage) LMM Formula:
#> LogRain ~ Gauge.Elevation + Steering.Wind.Speed + Total.Totals + 
#>     PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + 
#>     (1 | TrialDay)
#> 
#> Data subset used: oman [ Gauge.Day.Type == "Upwind"  &  Rain.Gauge.Measurement > 0 , ]
#> Number of observations: 1545, Number of groups: 292
#> 
#> Random effects:
#>  Groups   Name        Std.Dev.
#>  TrialDay (Intercept) 0.64489 
#>  Residual             1.26504 
#> 
#> Fixed effects:
#>               (Intercept)           Gauge.Elevation       Steering.Wind.Speed 
#>                   -1.4397                    0.4340                   -0.0964 
#>              Total.Totals       PC2.Dry.Temperature     PC1.Relative.Humidity 
#>                    0.0327                    0.1448                    0.1779 
#> PC1.Ground.Level.Pressure 
#>                   -0.0523 
#> 
#> ======================================================================
#> 
#> Downwind (Second Stage) LMM Formula:
#> LogRain - natural_pred ~ Gauge.Elevation + Target.H.01 + Target.H.02 + 
#>     Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + 
#>     Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + 
#>     Gauge.Elevation:Target.H.02 + (1 | TrialDay)
#> 
#> Data subset used: oman [ Gauge.Day.Type %in% c("Target", "Control")  &  Rain.Gauge.Measurement > 0 , ]
#> Number of observations: 4168, Number of groups: 488
#> 
#> Random effects:
#>  Groups   Name        Std.Dev.
#>  TrialDay (Intercept) 0.54789 
#>  Residual             1.35989 
#> 
#> Fixed effects:
#>                 (Intercept)             Gauge.Elevation 
#>                      0.3078                     -0.1958 
#>                 Target.H.01                 Target.H.02 
#>                      0.3157                      0.2398 
#>                 Target.H.03                 Target.H.04 
#>                      0.2368                     -0.1401 
#>                 Target.H.05                 Target.H.06 
#>                      0.4180                     -0.2002 
#>                 Target.H.07                 Target.H.08 
#>                      0.2309                      0.0786 
#>                 Target.H.09                 Target.H.10 
#>                      0.5648                      0.0527 
#> Gauge.Elevation:Target.H.01 Gauge.Elevation:Target.H.02 
#>                     -0.1335                     -0.1978 
summary(result)
#> Summary of Two Stage LMM Rainfall Enhancement Analysis
#> ======================================================================
#> 
#> Attribution Results (Assuming Log-Rainfall being Modelled):
#>     Estimate 95% Bootstrap CI Bootstrap P-Val Permutation P-Val
#> apo    6.27%               NA              NA                NA
#> apl    6.69%               NA              NA                NA
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> SATE Results:
#>             Estimate 95% Bootstrap CI Bootstrap P-Val Permutation P-Val
#> sate.mb       0.1153               NA              NA                NA
#> sate.ipw      0.0751               NA              NA                NA
#> sate.ipw.l    0.1132               NA              NA                NA
#> sate.ipw.ma   0.0813               NA              NA                NA
#> sate.aipw     0.0771               NA              NA                NA
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> 
#> ======================================================================
#> 
#> Upwind (First Stage) LMM:
#> Formula:
#> LogRain ~ Gauge.Elevation + Steering.Wind.Speed + Total.Totals + 
#>     PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + 
#>     (1 | TrialDay)
#> 
#> Data subset used: oman [ Gauge.Day.Type == "Upwind"  &  Rain.Gauge.Measurement > 0 , ]
#> Number of observations: 1545, Number of groups: 292
#> 
#> Random effects:
#>    Groups        Name  Variance
#>  TrialDay (Intercept) 0.4158771
#>  Residual             1.6003327
#> 
#> Fixed effects:
#>                           Estimate Std. Error t value
#> (Intercept)                -1.4397     0.3988 -3.6104
#> Gauge.Elevation             0.4340     0.0740  5.8631
#> Steering.Wind.Speed        -0.0964     0.0204 -4.7277
#> Total.Totals                0.0327     0.0083  3.9204
#> PC2.Dry.Temperature         0.1448     0.0527  2.7454
#> PC1.Relative.Humidity       0.1779     0.0223  7.9681
#> PC1.Ground.Level.Pressure  -0.0523     0.0201 -2.6041
#> 
#> ======================================================================
#> 
#> Downwind (Second Stage) LMM:
#> Formula:
#> LogRain - natural_pred ~ Gauge.Elevation + Target.H.01 + Target.H.02 + 
#>     Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + 
#>     Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + 
#>     Gauge.Elevation:Target.H.02 + (1 | TrialDay)
#> 
#> Data subset used: oman [ Gauge.Day.Type %in% c("Target", "Control")  &  Rain.Gauge.Measurement > 0 , ]
#> Number of observations: 4168, Number of groups: 488
#> 
#> Random effects:
#>    Groups        Name  Variance
#>  TrialDay (Intercept) 0.3001795
#>  Residual             1.8492958
#> 
#> Fixed effects:
#>                             Estimate Std. Error t value
#> (Intercept)                   0.3078     0.0626  4.9154
#> Gauge.Elevation              -0.1958     0.0606 -3.2293
#> Target.H.01                   0.3157     0.1362  2.3172
#> Target.H.02                   0.2398     0.1237  1.9389
#> Target.H.03                   0.2368     0.0925  2.5604
#> Target.H.04                  -0.1401     0.0890 -1.5738
#> Target.H.05                   0.4180     0.1319  3.1685
#> Target.H.06                  -0.2002     0.1493 -1.3413
#> Target.H.07                   0.2309     0.1876  1.2312
#> Target.H.08                   0.0786     0.1276  0.6163
#> Target.H.09                   0.5648     0.3062  1.8448
#> Target.H.10                   0.0527     0.1677  0.3141
#> Gauge.Elevation:Target.H.01  -0.1335     0.1437 -0.9287
#> Gauge.Elevation:Target.H.02  -0.1978     0.1196 -1.6539

# Extract quantities from the fitted downwind LMM
coef(result)$fixef_coef
#>                 (Intercept)             Gauge.Elevation 
#>                  0.30782064                 -0.19576293 
#>                 Target.H.01                 Target.H.02 
#>                  0.31568597                  0.23977373 
#>                 Target.H.03                 Target.H.04 
#>                  0.23675150                 -0.14009800 
#>                 Target.H.05                 Target.H.06 
#>                  0.41804072                 -0.20023955 
#>                 Target.H.07                 Target.H.08 
#>                  0.23091236                  0.07861614 
#>                 Target.H.09                 Target.H.10 
#>                  0.56479275                  0.05267170 
#> Gauge.Elevation:Target.H.01 Gauge.Elevation:Target.H.02 
#>                 -0.13348185                 -0.19775235 
head(coef(result)$ranef_coef)
#>          (Intercept)
#> 2013135 -0.191334499
#> 2013149  0.231232226
#> 2013150  0.841928530
#> 2013151  0.382101947
#> 2013152 -0.221780948
#> 2013157 -0.001806244
head(residuals(result))
#>          28        1690        1691        1709        1815        1833 
#> -0.86679304 -0.01061205  1.39145159 -0.33329984 -0.20167918  1.54411995 
head(fitted(result))
#>        28      1690      1691      1709      1815      1833 
#> 0.2336375 0.4155265 0.4233570 0.3881196 0.9763032 1.0005778 
varcomp(result)
#>                       TrialDay Residual
#> upwind_lmm           0.4158771 1.600333
#> downwind_lmm         0.3001795 1.849296
#> downwind_target_lmm  0.3821506 1.811551
#> downwind_control_lmm 0.3577615 1.781349
head(predict(result))
#>        28      1690      1691      1709      1815      1833 
#> 0.2336375 0.4155265 0.4233570 0.3881196 0.9763032 1.0005778 

# Model diagnostic plots for the downwind LMM
plot(result, plot_type = "model")
#> Registered S3 method overwritten by 'car':
#>   method           from
#>   na.action.merMod lme4

```
