# Focused validation tests for rain_attr()
#
make_validation_rain_attr <- function(
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
    attr_type = "ThoEtAl",
    x_downwind_name = "Gauge.Elevation",
    bootstrap = FALSE,
    bootstrap_option = bootstrap_opt()) {

  rain_attr(
    data = oman,
    upwind_lmm_formula =
      LogRain ~ Gauge.Elevation + Steering.Wind.Speed + Total.Totals +
      PC2.Dry.Temperature + PC1.Relative.Humidity +
      PC1.Ground.Level.Pressure + (1 | TrialDay),
    instr_pred_name = instr_pred_name,
    instr_pred_type = instr_pred_type,
    downwind_lmm_formula = downwind_lmm_formula,
    downwind_logistic_formula = downwind_logistic_formula,
    downwind_propensity_formula = downwind_propensity_formula,
    rain_col_name = "Rain.Gauge.Measurement",
    upwind_subset = Gauge.Day.Type == "Upwind",
    downwind_subset = Gauge.Day.Type %in% c("Target", "Control"),
    downwind_target_subset = Gauge.Day.Type == "Target",
    downwind_control_subset = Gauge.Day.Type == "Control",
    positive_subset = Rain.Gauge.Measurement > 0,
    attr_type = attr_type,
    x_downwind_name = x_downwind_name,
    bootstrap = bootstrap,
    bootstrap_option = bootstrap_option,
    permutation = FALSE
  )
}


test_that("rain_attr validates instrumental prediction specification", {
  expect_error(
    make_validation_rain_attr(
      instr_pred_name = "not_in_downwind_formula"
    )
  )

  expect_error(
    make_validation_rain_attr(
      instr_pred_type = "invalid"
    )
  )
})


test_that("rain_attr validates treatment and positive-rainfall definitions", {
  expect_error(
    make_validation_rain_attr(
      downwind_propensity_formula =
        (Gauge.Day.Type == "Control") ~ Total.Totals +
        PC1.Dry.Temperature + PC1.Relative.Humidity +
        PC1.Ground.Level.Pressure
    )
  )

  expect_error(
    make_validation_rain_attr(
      downwind_logistic_formula =
        (Rain.Gauge.Measurement >= 0) ~ Gauge.Elevation + natural_pred
    )
  )
})


test_that("rain_attr validates attribution specification", {
  expect_error(
    make_validation_rain_attr(
      attr_type = "invalid"
    )
  )

  expect_error(
    make_validation_rain_attr(
      x_downwind_name = "not_a_variable"
    )
  )
})


test_that("rain_attr validates key bootstrap requirements", {
  expect_error(
    make_validation_rain_attr(
      bootstrap = TRUE,
      bootstrap_option = bootstrap_opt(
        B_bootstrap = 1,
        bootstrap_zero = TRUE
      )
    )
  )

  expect_error(
    make_validation_rain_attr(
      bootstrap = TRUE,
      bootstrap_option = bootstrap_opt(
        bootstrap_zero = FALSE,
        CI_level = 1.1
      )
    )
  )

  expect_error(
    make_validation_rain_attr(
      bootstrap = TRUE,
      bootstrap_option = bootstrap_opt(
        bootstrap_zero = FALSE,
        bootstrap_type = 'invalid_type'
      )
    )
  )
})
