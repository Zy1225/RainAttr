# Daily Oman Ionizer Operation Schedule

This dataset contains the daily operation schedule of the ionizers
during the 2013–2018 Oman rainfall enhancement trial, as used in the
analyses of Chambers et al. (2022a,b).

## Usage

``` r
ionizer_operation
```

## Format

A data frame with 740 rows (days) and 12 variables:

- TrialDay:

  Identifier for each day in the six-year trial

- Year:

  Character variable indicating the year

- H1, H2, H3, H4, H5, H6, H7, H8, H9, H10:

  Binary indicator variables for the operating status of the ionizers: 1
  = ionizer active on that day, 0 = ionizer inactive, NA = ionizer not
  yet deployed

## References

- Chambers, R., Beare, S., Peak, S. and Al-Kalbani, M. (2022a) Nudging a
  Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall
  Enhancement Trial in Oman. *International Statistical Review*, 90:
  346–373.

- Chambers, R., Ranjbar, S., Salvati, N., and Pacini, B. (2022b)
  Weighting, Informativeness and Causal Inference, with an Application
  to Rainfall Enhancement. *Journal of the Royal Statistical Society
  Series A: Statistics in Society*, 185: 1584–1612
