# Package index

## Main function

Main function for fitting the two-stage model and performing estimation
and inference of attribution and sample average treatment effect.

- [`rain_attr()`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
  : Attribution and Sample Average Treatment Effect for Rainfall
  Enhancement Trial Data

## Fitted objects and methods

Documentation for fitted rain_attr objects and their associated S3
methods, including coefficients, residuals, predictions, variance
components, diagnostic plots, and summaries.

- [`coef(`*`<rain_attr>`*`)`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md)
  [`print(`*`<rain_attr>`*`)`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md)
  [`residuals(`*`<rain_attr>`*`)`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md)
  [`fitted(`*`<rain_attr>`*`)`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md)
  [`varcomp(`*`<rain_attr>`*`)`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md)
  [`predict(`*`<rain_attr>`*`)`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md)
  [`plot(`*`<rain_attr>`*`)`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md)
  [`summary(`*`<rain_attr>`*`)`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md)
  [`print(`*`<summary.rain_attr>`*`)`](https://zy1225.github.io/RainAttr/reference/rain_attr-class.md)
  : Class "rain_attr" of Two-Stage LMM Fitted to Rainfall Enhancement
  Trial Data

## Inference options

Functions for specifying bootstrap and permutation settings.

- [`bootstrap_opt()`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md)
  : Bootstrap Options
- [`permutation_opt()`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md)
  : Permutation Options

## Inference procedures

Detailed documentation of the bootstrap and permutation procedures used
for inference in rainfall enhancement trials.

- [`bootstrap_downwind()`](https://zy1225.github.io/RainAttr/reference/bootstrap_downwind.md)
  : Bootstrap Procedure for Rainfall Enhancement Trial Data
- [`permutation_ionizer()`](https://zy1225.github.io/RainAttr/reference/permutation_ionizer.md)
  : Permutation-Based Procedure for Rainfall Enhancement Trial Data

## Exploratory analysis

- [`eda()`](https://zy1225.github.io/RainAttr/reference/eda.md) :
  Exploratory Data Analysis for Rainfall Enhancement Trial Data

## Example data

- [`oman`](https://zy1225.github.io/RainAttr/reference/oman.md) : Oman
  Rainfall Enhancement Trial Dataset
- [`gaugeday_downwind`](https://zy1225.github.io/RainAttr/reference/gaugeday_downwind.md)
  : Gauge-Day Downwind Indicator
- [`ionizer_location`](https://zy1225.github.io/RainAttr/reference/ionizer_location.md)
  : Oman Ionizer Location and Deployment Timing
- [`ionizer_operation`](https://zy1225.github.io/RainAttr/reference/ionizer_operation.md)
  : Daily Oman Ionizer Operation Schedule
