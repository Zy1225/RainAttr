# Oman Rainfall Enhancement Trial Dataset

This dataset contains gauge-day observations from the 2013–2018 Oman
rainfall enhancement trial, which was used in the analyses of Chambers
et al. (2022a,b).

## Usage

``` r
oman
```

## Format

A data frame with 122259 rows (gauge-day observations) and 55 variables:

- Year, YearDay, Month:

  Time-related variables

- TrialDay:

  Identifier for each day in the six-year trial

- Year...2013, Year...2014, Year...2015, Year...2016, Year...2017,
  Year...2018:

  Binary indicator variables for each year in the trial period

- H1on, H2on, H3on, H4on, H5on, H6on, H7on, H8on, H9on, H10on:

  Binary indicator variables for the operating status of ionizers (1 =
  on, 0 = off), with NA indicating ionizers not yet deployed

- Gauge.ID, Gauge.Latitude, Gauge.Longitude, Gauge.Elevation:

  Gauge identifiers, location and elevation

- Elevated.Gauge:

  Binary indicator variable for gauges located above 1km elevation (1 =
  elevated, 0 = not elevated)

- Gauge.Elevation...1km:

  Gauge elevation for non-elevated gauges, computed as Gauge.Elevation ×
  (1 - Elevated.Gauge)

- Gauge.Elevation...1km.1:

  Gauge elevation for elevated gauges, computed as Gauge.Elevation ×
  Elevated.Gauge

- Rain.Gauge.Measurement:

  Observed rainfall at the raw scale

- Rainfall.Event:

  Binary indicator variable for a rainfall event (1 = rainfall, 0 = no
  rainfall)

- Positive.Rainfall:

  Observed rainfall for rainfall events; NA if no rainfall occurred

- LogRain:

  Log-transformed rainfall for rainfall events; NA if no rainfall
  occurred

- Rainfall.Measurement.Status:

  Character variable indicating the gauge’s location relative to
  deployed ionizers based on the day’s wind direction, with possible
  values: "Upwind", "Downwind", "Out of Scope"

- Target.H.01, Target.H.02, Target.H.03, Target.H.04, Target.H.05,
  Target.H.06, Target.H.07, Target.H.08, Target.H.09, Target.H.10:

  Binary indicator variables showing whether the gauge is downwind of
  the ionizer and the ionizer has been deployed and turned on on that
  day (1 = yes, 0 = no)

- Gauge.Day.Type:

  Character variable that further classifies "Downwind" observations -
  "Target" – downwind of at least one active ionizer (i.e., at least one
  of Target.H.01 – Target.H.10 = 1) - "Control" – downwind of all
  inactive ionizers (i.e., all Target.H.01 – Target.H.10 = 0)

  Observations that are "Upwind" or "Out of Scope" retain the same
  classification as Rainfall.Measurement.Status. Thus, the variable has
  four possible values: "Upwind", "Target", "Control", "Out of Scope".

- PC1.Dry.Temperature, PC2.Dry.Temperature,PC1.Relative.Humidity,
  PC2..Relative.Humidity,PC1.Ground.Level.Pressure:

  Principal component of meteorological variables

- Steering.Wind.Direction, Steering.Wind.Principal.Direction,
  Steering.Wind.Speed:

  Wind-related measures

- Lifted.Index, Total.Totals:

  Atmospheric instability indices

- LCL.Pressure:

  A measure of cloud base height

- Precipitable.Water:

  Amount of rainfall if a column of the atmosphere were to be
  precipitated

## References

- Chambers, R., Beare, S., Peak, S. and Al-Kalbani, M. (2022a) Nudging a
  Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall
  Enhancement Trial in Oman. *International Statistical Review*, 90:
  346–373.

- Chambers, R., Ranjbar, S., Salvati, N., and Pacini, B. (2022b)
  Weighting, Informativeness and Causal Inference, with an Application
  to Rainfall Enhancement. *Journal of the Royal Statistical Society
  Series A: Statistics in Society*, 185: 1584–1612
