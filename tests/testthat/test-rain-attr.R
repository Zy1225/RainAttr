fit_oman <- function(attr_type) {
  rain_attr(
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
    attr_type = attr_type,
    x_downwind_name = "Gauge.Elevation",
    bootstrap = FALSE,
    permutation = FALSE
  )
}

expected_sate <- c(
  sate.mb = 0.11531,
  sate.ipw = 0.07513,
  sate.ipw.l = 0.11324,
  sate.ipw.ma = 0.08135,
  sate.aipw = 0.07706
)


test_that("rain_attr reproduces standard Oman results for attr_type = ThoEtAl", {
  result <- fit_oman("ThoEtAl")

  expect_s3_class(result, "rain_attr")
  expect_equal(result$hatattr$apo, 0.06267536, tolerance = 1e-6)
  expect_equal(result$hatattr$apl, 0.06686623, tolerance = 1e-6)
  expect_equal(unlist(result$hatsate), expected_sate, tolerance = 1e-3)
})

test_that("rain_attr reproduces standard Oman results for attr_type = ChambersEtAl attr_type", {
  result <- fit_oman("ChambersEtAl")

  expect_equal(result$hatattr$apo, 0.1018356, tolerance = 1e-6)
  expect_equal(result$hatattr$apl, 0.1133819, tolerance = 1e-6)
  expect_equal(unlist(result$hatsate), expected_sate, tolerance = 1e-3)
})

test_that("rain_attr reproduces standard Oman results for attr_type = ChambersEtAl_No_Winsorize", {
  result <- fit_oman("ChambersEtAl_No_Winsorize")

  expect_equal(result$hatattr$apo, 0.1018484, tolerance = 1e-6)
  expect_equal(result$hatattr$apl, 0.1133978, tolerance = 1e-6)
  expect_equal(unlist(result$hatsate), expected_sate, tolerance = 1e-3)
})

test_that("rain_attr reproduces standard Oman results for attr_type = No ", {
  result <- fit_oman("No")

  expect_equal(result$hatattr$apo, 0.05953303, tolerance = 1e-6)
  expect_equal(result$hatattr$apl, 0.06330156, tolerance = 1e-6)
  expect_equal(unlist(result$hatsate), expected_sate, tolerance = 1e-3)
})
