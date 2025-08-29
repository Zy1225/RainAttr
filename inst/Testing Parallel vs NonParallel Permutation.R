rm(list=ls())
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
                 attr_type = 'Proposed',
                 x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
                 target_only = FALSE,
                 permutation = TRUE,
                 permutation_option = permutation_option(
                   B_permutation = 500,
                   permute_between_ionizer = T,
                   permute_all_ionizers_between_day = T,
                   permute_between_gaugeday = T,
                   ionizer_operation_input = ionizer_operation,
                   gaugeday_downwind_input = gaugeday_downwind,
                   year_ionizer_list =
                     list(
                       '2013' = c('H1','H2'),
                       '2014' = c('H1','H2','H3','H4'),
                       '2015' = c('H1','H2','H3','H4','H5','H6'),
                       '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
                       '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
                       '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
                     ),
                   data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
                   ionizer_operation_year_column_name = 'Year',
                   ionizer_operation_day_column_name = 'TrialDay',
                   permutation_seed = 123,
                   permutation_parallel = F,
                   permutation_parallel_num_worker = NULL
                 ))

end.time = Sys.time()

colMeans(asd$permutation_result$hatattr)
# apo        apl
# 0.01833859 0.01916830
colMeans(asd$permutation_result$hatsate)
# sate.mb    sate.ipw  sate.ipw.l sate.ipw.ma   sate.aipw
# 0.02785166  0.01697089  0.02775810  0.01154884  0.01678482

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
                          attr_type = 'Proposed',
                          x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
                          target_only = FALSE,
                          permutation = TRUE,
                          permutation_option = permutation_option(
                            B_permutation = 500,
                            permute_between_ionizer = T,
                            permute_all_ionizers_between_day = T,
                            permute_between_gaugeday = T,
                            ionizer_operation_input = ionizer_operation,
                            gaugeday_downwind_input = gaugeday_downwind,
                            year_ionizer_list =
                              list(
                                '2013' = c('H1','H2'),
                                '2014' = c('H1','H2','H3','H4'),
                                '2015' = c('H1','H2','H3','H4','H5','H6'),
                                '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
                                '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
                                '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
                              ),
                            data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
                            ionizer_operation_year_column_name = 'Year',
                            ionizer_operation_day_column_name = 'TrialDay',
                            permutation_seed = 123,
                            permutation_parallel = T,
                            permutation_parallel_num_worker = parallel::detectCores() - 1
                          ))
end.parallel.time = Sys.time()

colMeans(asd_parallel$permutation_result$hatattr)
# apo        apl
# 0.01832737 0.01911959
colMeans(asd_parallel$permutation_result$hatsate)
# sate.mb    sate.ipw  sate.ipw.l sate.ipw.ma   sate.aipw
# 0.02798632  0.01790311  0.02789781  0.01287128  0.01759901




#
asd$hatattr$apo; asd$hatattr$apl
asd_parallel$hatattr$apo; asd_parallel$hatattr$apl

unlist(asd$hatsate)
unlist(asd_parallel$hatsate)

end.time - start.time
end.parallel.time - start.parallel.time
