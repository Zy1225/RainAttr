test_that("eda partitions downwind observations completely into target and control", {
  x <- eda(
    eda_type = "num_obs_days",
    data = oman,
    day_column_name = "TrialDay",
    upwind_subset = Gauge.Day.Type == "Upwind",
    downwind_subset = Gauge.Day.Type %in% c("Target", "Control"),
    downwind_target_subset = Gauge.Day.Type == "Target",
    downwind_control_subset = Gauge.Day.Type == "Control",
    positive_subset = Rain.Gauge.Measurement > 0
  )

  expect_equal(
    unname(x$num_obs[2, ]),
    unname(x$num_obs[3, ] + x$num_obs[4, ])
  )
})
