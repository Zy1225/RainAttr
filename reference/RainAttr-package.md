# RainAttr: Attribution and Sample Average Treatment Effect for Rainfall Enhancement Trial Data

`RainAttr` provides tools for estimation and inference of attribution
and sample average treatment effect in rainfall enhancement trials using
a two-stage linear mixed model approach based on Chambers et al. (2022).

## Details

The main functions are:

- [`rain_attr`](https://zy1225.github.io/RainAttr/reference/rain_attr.md)
  for fitting the two-stage model, estimating attribution and SATEs, and
  optionally performing bootstrap-based and permutation-based inference.

- [`bootstrap_opt`](https://zy1225.github.io/RainAttr/reference/bootstrap_opt.md)
  for configuring bootstrap inference.

- [`permutation_opt`](https://zy1225.github.io/RainAttr/reference/permutation_opt.md)
  for configuring permutation inference.

- [`eda`](https://zy1225.github.io/RainAttr/reference/eda.md) for
  exploratory data analysis.

See the package vignette for a complete worked example based on the Oman
rainfall enhancement trial.

## References

- Chambers, R., Beare, S., Peak, S. and Al-Kalbani, M. (2022) Nudging a
  Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall
  Enhancement Trial in Oman. *International Statistical Review*, 90:
  346–373.

## See also

Useful links:

- <https://github.com/Zy1225/RainAttr>

- <https://zy1225.github.io/RainAttr/>

- Report bugs at <https://github.com/Zy1225/RainAttr/issues>

## Author

**Maintainer**: Zhi Yang Tho <zhiyang.tho@anu.edu.au>

Authors:

- Raymond Chambers <Raymond.Chambers@anu.edu.au>

- A. H. Welsh <Alan.Welsh@anu.edu.au>
