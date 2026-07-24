# Class "rain_attr" of Two-Stage LMM Fitted to Rainfall Enhancement Trial Data

Objects of class `rain_attr` represent results from the two-stage
LMM-based rainfall enhancement analysis, created by calls to
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).

## Usage

``` r
# S3 method for class 'rain_attr'
coef(object, model = "downwind_lmm", ...)

# S3 method for class 'rain_attr'
print(object, ...)

# S3 method for class 'rain_attr'
residuals(
  object,
  model = "downwind_lmm",
  residual_type = NULL,
  residual_scaled = TRUE,
  ...
)

# S3 method for class 'rain_attr'
fitted(object, model = "downwind_lmm", ...)

# S3 method for class 'rain_attr'
varcomp(object, ...)

# S3 method for class 'rain_attr'
predict(
  object,
  newdata = NULL,
  model = "downwind_lmm",
  re_include = TRUE,
  fixef_include = TRUE,
  allow.new.levels = FALSE,
  predict_type = "link",
  ...
)

# S3 method for class 'rain_attr'
plot(
  object,
  plot_type = c("bootstrap", "permutation"),
  plot_quantity = c("attr", "sate"),
  model = "downwind_lmm",
  residual_type = NULL,
  residual_scaled = TRUE,
  re_include = TRUE,
  fixef_include = TRUE,
  allow.new.levels = FALSE,
  predict_type = "link",
  ...
)

# S3 method for class 'rain_attr'
summary(object, ...)

# S3 method for class 'summary.rain_attr'
print(summary_object, ...)
```

## Arguments

- object:

  An **R** object of class `rain_attr`, i.e., an output from
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).

- model:

  An optional character string specifying which fitted model to focus
  on. Must be one of `"upwind_lmm"`, `"downwind_lmm"` (default),
  `"downwind_target_lmm"`, `"downwind_control_lmm"`,
  `"downwind_logistic"`, or `"downwind_propensity"`. See "Details" for
  more information.

- ...:

  Additional arguments passed to methods.

- residual_type:

  An optional character string specifying the type of residuals. Must be
  one of `"response"` (default for LMMs), `"deviance"` (default for
  GLMs), `"working"`, or `"pearson"`. An additional choice for GLMs is
  `"partial"`, but this residual type can only be used with
  `residuals.rain_attr` but not `plot.rain_attr`. See
  [`residuals.merMod`](https://rdrr.io/pkg/lme4/man/residuals.merMod.html)
  and [`residuals.glm`](https://rdrr.io/r/stats/glm.summaries.html) for
  definition of each residual type of LMMs and GLMs, respectively.

- residual_scaled:

  An optional logical. Only used when `model` is one of the LMMs. If
  `TRUE` (default), the residuals are scaled by residual standard
  deviation i.e., \\\hat{\sigma}\_e\\. If `FALSE`, no scaling is done
  for the residuals.

- newdata:

  An optional data frame in which to look for variables to predict from,
  based on the chosen `model`. If omitted, predictions are based on the
  original data used to fit `model`.

- re_include:

  An optional logical. Only used when `model` is one of the LMMs. If
  `TRUE` (default), random intercepts are included when predicting. If
  `FALSE`, random intercepts are not included when predicting.

- fixef_include:

  An optional logical. Only used when `model` is one of the LMMs. If
  `TRUE` (default), fixed effects are included in predicting. If
  `FALSE`, fixed effects are not included in predicting.

- allow.new.levels:

  An optional logical. Only used when `model` is one of the LMMs. If
  `TRUE`, new days (groups) are allowed in `newdata`, where the
  prediction will use zeros (population-level mean) as the predicted
  random effects. If `FALSE` (default), such new days (groups) in
  `newdata` will trigger an error.

- predict_type:

  An optional character string specifying the type of prediction. Only
  used when `model` is one of the GLMs. Must be one of `"link"`
  (default), `"response"`, or `"terms"`. `"link"` gives predictions on
  the scale of the linear predictors of the GLMs i.e., log-odds.
  `"response"` gives predictions on the scale of the response variable,
  i.e., the predicted probabilities. `"terms"` return a matrix giving
  the fitted values of each term in the model formula on the linear
  predictor (log-odds) scale, so `"terms"` can only be used with
  `predict.rain_attr` but not `plot.rain_attr`.

- plot_type:

  An optional character string or vector specifying the type of plots.
  Can be a single value: `"bootstrap"`, `"permutation"`, or `"model"`,
  or a vector of length 2: `c("bootstrap", "permutation")` (default).
  `"bootstrap"` or `"permutation"` should only be used when the original
  call to
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
  in creating `object` has `bootstrap = TRUE` or `permutation = TRUE`,
  respectively, which will plot the bootstrap or permutation
  distribution of the chosen `plot_quantity`. `"model"` provides two
  diagnostic plots for the chosen `model`: a plot of residuals against
  fitted values, and a Q-Q plot of residuals.

- plot_quantity:

  An optional character string or vector specifying the quantity whose
  bootstrap or permutation distribution is to be plotted. Only used when
  `plot_type` includes `"bootstrap"` or `"permutation"`. Can be a single
  value: `"attr"` or `"sate"`, or a vector of length 2:
  `c("attr","sate")` (default). `"attr"` stands for the attribution
  estimates (`apo`, `apl`), while `"sate"` stands for the sample average
  treatment effect estimates (`sate.mb`, `sate.ipw`, `sate.ipw.l`,
  `sate.ipw.ma`, `sate.aipw`). See
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
  for more details.

- summary_object:

  An **R** object of class `summary.rain_attr`, i.e., an output from
  applying `summary.rain_attr` to an **R** object of class `rain_attr`.

## Details

A `rain_attr` object contains fitted LMMs:

- `upwind_lmm`:

  Upwind (first stage) LMM.

- `downwind_lmm`:

  Downwind (second stage) LMM.

- `downwind_target_lmm`:

  Downwind (second stage) treatment-only LMM.

- `downwind_control_lmm`:

  Downwind (second stage) control-only LMM.

and fitted GLMs:

- `downwind_logistic`:

  Downwind (second stage) logistic model of rainfall event indicator.

- `downwind_propensity`:

  Downwind (second stage) propensity score model for the treatment
  indicator.

The object also includes attribution estimates, SATE, results from
bootstrap and permutation-based procedures if `bootstrap = TRUE` or
`permutation = TRUE` were specified in the original
[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
call, the data supplied in the original call with an additional column
containing the fitted values generated from the upwind (first stage)
LMM., and the specified arguments in the original call.

## S3 Methods

- `coef`:

  Computes the fixed effect coefficients and the empirical best linear
  unbiased predictors (EBLUPs) of random intercepts if the selected
  `model` is one of the LMMs, or the regression coefficients if the
  selected `model` is one of the GLMs.

- `print`:

  Print the point estimates of attributions and SATE and indicate
  whether bootstrap and permutation has been performed. Also print the
  model formula, data subset, number of observations, number of days
  (groups), and point estimates of fixed effect coefficients for the
  upwind (first stage) and downwind (second stage) LMMs.

- `residuals`:

  Compute residuals of the selected `model`. Residual type is controlled
  by `residual_type`, while scaling (only for LMMs) is controlled by
  `residual_scaled`.

- `fitted`:

  Compute fitted values of the selected `model`. For more flexible
  access to fitted values, use `predict.rain_attr`.

- `varcomp`:

  Compute the random intercept variance and error variance estimates for
  each of the four LMMs.

- `predict`:

  Compute predictions from the selected `model`. The `newdata` argument
  allows predictions for a new data frame; if omitted, predictions are
  based on the original data used to fit the model. For GLMs, the type
  of prediction is controlled by `predict_type`. For LMMs, `re_include`,
  `fixef_include`, `allow.new.levels` control how fixed effects and
  random effects are handled when making predictions.

- `plot`:

  Depending on `plot_type`, it can plot bootstrap and/or permutation
  distributions of attribution and/or SATE estimates, or model
  diagnostic plots (residuals vs fitted and Q-Q plots) for a selected
  `model`. The `plot_quantity` argument controls which estimates to
  display for bootstrap and/or permutation plots. For bootstrap plots, a
  solid vertical line represents the estimate based on the original data
  and a dotted vertical line represents zero. For permutation plots, a
  solid vertical line represents the estimate based on the original
  data. For diagnostic plots of `model = "downwind_lmm"`, the treatment
  and control observations are plotted using different colors and
  shapes. Diagnostic plots for GLM-type `model` should be interpreted
  with caution. In particular, logistic regressions (for both
  `"downwind_logistic"` and `"downwind_propensity"`) are inherently
  curvilinear, which often leads to unusual patterns in the
  residuals–fitted values plot even when the model is correctly
  specified. Moreover, the validity of GLMs does not require normally
  distributed residuals.

- `summary`:

  Returns a list of summary statistics of the two-stage LMM approach,
  which is an **R** object of class `summary.rain_attr`. For details of
  this list, see the "Value from `summary.rain_attr`" section below.

- `print.summary`:

  Print the point estimates of attributions and SATE, along with their
  bootstrap confidence intervals, bootstrap p-values and permutation
  p-values if bootstrap and permutation has been performed. Also print
  the model formula, data subset, number of observations, number of
  unique days (groups), variance component estimates, and fixed effect
  coefficient estimates along with their standard error and t-values
  computed from
  [`summary.merMod`](https://rdrr.io/pkg/lme4/man/summary.merMod.html).

## Value from `summary.rain_attr`

The following describes the output from `summary.rain_attr`:

- `attr_table`:

  A data frame containing attribution estimates, along with their
  bootstrap confidence intervals, bootstrap p-values, and permutation
  p-values. If no bootstrap or permutation has been carried out, their
  corresponding confidence intervals or p-values are set as `NA`.

- `sate_table`:

  A data frame containing SATE estimates, along with their bootstrap
  confidence intervals, bootstrap p-values, and permutation p-values. If
  no bootstrap or permutation has been carried out, their corresponding
  confidence intervals or p-values are set as `NA`.

- `data_name`:

  The name of the data frame used in the original call to
  [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md).

- `upwind_subset_expr`, `downwind_subset_expr`:

  Expressions used to subset the data for fitting the upwind and
  downwind LMMs.

- `upwind_formula`, `downwind_formula`:

  Model formulas used for fitting the upwind (first stage) and downwind
  (second stage) LMMs.

- `upwind_n_obs`, `downwind_n_obs`:

  Number of observations used for fitting the upwind and downwind LMMs.

- `upwind_n_groups`, `downwind_n_groups`:

  Number of unique days (groups) used for fitting the upwind and
  downwind LMMs.

- `upwind_summary`, `downwind_summary`:

  Objects of class
  [`summary.merMod`](https://rdrr.io/pkg/lme4/man/merMod-class.html)
  containing summaries of the fitted upwind and downwind LMMs.

- `upwind_fitted`, `downwind_fitted`:

  Fitted values of the upwind and downwind LMMs, which include both the
  fixed effects and the predicted random intercepts.

- `upwind_residuals`, `downwind_residuals`:

  Scaled residuals of the upwind and downwind LMMs, which take the
  general form of (observed - fitted)/\\\hat{\sigma}\_e\\, where the
  fitted values are either `upwind_fitted` or `downwind_fitted` and
  \\\hat{\sigma}\_e\\ denote the corresponding estimated error variance.

- `upwind_lmm_fixef`, `downwind_lmm_fixef`:

  Matrices of fixed effect coefficients for the upwind and downwind
  LMMs, along with their corresponding standard error estimates and
  t-value computed from
  [`summary.merMod`](https://rdrr.io/pkg/lme4/man/summary.merMod.html).

- `upwind_lmm_varcomp`, `downwind_lmm_varcomp`:

  Data frames containing the estimated variance components for the
  upwind and downwind LMMs.

## See also

[`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
