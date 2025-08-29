rm(list=ls()); gc()
load('data/oman.rda')
load('data/ionizer_operation.rda')
load('data/gaugeday_downwind.rda')
source('R/Functions.R')

start.time = Sys.time()
asd  = rain_attr(data = oman,
                 upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
                 instr_pred_name = 'natural_pred',
                 instr_pred_type = 'Unconditional',
                 downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
                 downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
                 downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
                 rain_col_name = 'Rain.Gauge.Measurement',
                 upwind_subset = Gauge.Day.Type == 'Upwind',
                 downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
                 downwind_target_subset = Gauge.Day.Type == 'Target',
                 downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
                 attr_type = 'ThoEtAl',
                 x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
                 target_only = FALSE,
                 bootstrap = TRUE,
                 bootstrap_option = bootstrap_option(
                   B_bootstrap = 10000,
                   bootstrap_seed = 123,
                   bootstrap_parallel = F,
                   bootstrap_parallel_num_worker = NULL
                 ))

end.time = Sys.time()

colMeans(asd$bootstrap_result$hatattr)
# apo        apl
# 0.06423675 0.06920307
colMeans(asd$bootstrap_result$hatsate)
# sate.mb    sate.ipw  sate.ipw.l sate.ipw.ma   sate.aipw
# 0.1038889   0.1151946   0.1036914   0.1048006   0.1182679

start.parallel.time = Sys.time()
asd_parallel  = rain_attr(data = oman,
                          upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
                          instr_pred_name = 'natural_pred',
                          instr_pred_type = 'Unconditional',
                          downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
                          downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
                          downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
                          rain_col_name = 'Rain.Gauge.Measurement',
                          upwind_subset = Gauge.Day.Type == 'Upwind',
                          downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
                          downwind_target_subset = Gauge.Day.Type == 'Target',
                          downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
                          attr_type = 'ThoEtAl',
                          x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
                          target_only = FALSE,
                          bootstrap = TRUE,
                          bootstrap_option = bootstrap_option(
                            B_bootstrap = 10000,
                            bootstrap_seed = 123,
                            bootstrap_parallel = T,
                            bootstrap_parallel_num_worker = parallel::detectCores() - 1
                          ))
end.parallel.time = Sys.time()

colMeans(asd_parallel$bootstrap_result$hatattr)
# apo        apl
# 0.06447918 0.06947697
colMeans(asd_parallel$bootstrap_result$hatsate)
# sate.mb    sate.ipw  sate.ipw.l sate.ipw.ma   sate.aipw
# 0.1044085   0.1159940   0.1042058   0.1056489   0.1190891




#
asd$hatattr$apo; asd$hatattr$apl
asd_parallel$hatattr$apo; asd_parallel$hatattr$apl

unlist(asd$hatsate)
unlist(asd_parallel$hatsate)

end.time - start.time #Time difference of 2.596566 hours
end.parallel.time - start.parallel.time #Time difference of 24.24193 mins

save.image('inst/result_parallel_vs_nonparallel_bootstrap.Rdata')
