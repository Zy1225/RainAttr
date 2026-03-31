# R script to demonstrate functionality of RainAttr package


rm(list=ls())
load('data/oman.rda')
load('data/ionizer_operation.rda')
load('data/gaugeday_downwind.rda')
load('data/ionizer_location.rda')
devtools::load_all(".")

#### To replicate point estimation results of "Weighting, informativeness and causal inference, with an application to rainfall enhancement" by Chambers et al. (2021) in JRSSA ####

replicate_result = rain_attr(
  #input data
  data = oman,

  #LMM formula fitted to the upwind observations
  upwind_lmm_formula = LogRain ~ Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),

  #Variable name for storing instrumental prediction from upwind LMM
  instr_pred_name = 'natural_pred',

  #Types of instrumental prediction: 'Unconditional' or 'Conditional'
  instr_pred_type = 'Unconditional',

  #LMM formula fitted to the downwind observations
  downwind_lmm_formula = LogRain ~ Year...2013 + Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation...1km:Target.H.01 + Gauge.Elevation...1km:Target.H.02 + + Gauge.Elevation...1km.1:Target.H.01 + Gauge.Elevation...1km.1:Target.H.02 + (1|TrialDay),

  #Formula for fitting a logistic regression model to the indicator of rainfall event (used for bootstrapping of rainfall events)
  downwind_logistic_formula = NULL,

  #Formula for fitting the propensity score model to downwind observations, where responses are indicators of whether each downwind observations is a 'Target'
  downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,

  #Specify the column in input data that contains the raw rainfall
  rain_col_name = 'Rain.Gauge.Measurement',

  #Logical expression identifying subset of observations to which upwind_lmm_formula is fitted
  upwind_subset = Gauge.Day.Type == 'Upwind' & Year!= 2013,

  #Logical expression identifying subset of observations to which downwind_lmm_formula is fitted
  downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),

  #Logical expression identifying subset of observations of those that are 'Target'
  downwind_target_subset = Gauge.Day.Type == 'Target',

  #Logical expression identifying subset of observations of those that are 'Control'
  downwind_control_subset = Gauge.Day.Type == 'Control',

  #Logical expression identifying subset of observations with rainfall event
  positive_subset = Rain.Gauge.Measurement > 0,

  #Types of correction (for back-transformation bias) used for computing the attribution estimate: 'Chambers_Chandra', 'ThoEtAl', or 'No'
  attr_type = 'ChambersEtAl',

  #Vector of variable names in downwind_lmm_formula to identify non-ionizer related covariates, whose effects are not included in the calculation of attribution and SATE
  x_downwind_name = c('Year...2013' , 'Year...2014' , 'Year...2016' , 'Year...2017' , 'Year...2018', 'Gauge.Elevation...1km', 'Gauge.Elevation...1km.1', 'natural_pred'),

  #Logical indicating whether the attribution estimates should be based on only 'Target' observations (TRUE), or based on both 'Target' and 'Control' observations (FALSE)
  target_only = FALSE,

  #Indicator for whether to perform bootstrap inference on attribution and SATE
  bootstrap = FALSE,

  #Specification of various option for bootstrap, e.g., bootstrap_opt(bootstrap_type = 'PREB1'), currently support bootstrap_type = 'PREB0', 'PREB1', 'PREB2', 'REB0', 'REB1', 'REB2', 'MREB1'
  bootstrap_option = bootstrap_opt(),

  #Indicator for whether to perform permutation inference on attribution and SATE
  permutation = FALSE,

  #Specification of various option for permutation,  e.g., whether to permute the operating states between ionizers, between days, between gaugedays
  permutation_option = permutation_opt()
  )

#Table 1:
replicate_result$all_fitted_models$downwind_propensity_fit

#Table 2:
lme4::fixef(replicate_result$all_fitted_models$downwind_lmm_fit)

#Table 3:
as.data.frame(lme4::VarCorr(replicate_result$all_fitted_models$downwind_lmm_fit))

#Table 4:
replicate_result$all_fitted_models$downwind_positive_target_lmm_fit
replicate_result$all_fitted_models$downwind_positive_control_lmm_fit

#Table 5
as.data.frame(lme4::VarCorr(replicate_result$all_fitted_models$downwind_positive_target_lmm_fit))
as.data.frame(lme4::VarCorr(replicate_result$all_fitted_models$downwind_positive_control_lmm_fit))

#Table 6: the row corresponding to Estimate for LogRain - currently do not support exponentiation of the SATE estimate
unlist(replicate_result$hatsate)

#Attribution estimates
c(replicate_result$hatattr$apl, replicate_result$hatattr$apo)

#### Documentations ####

#Two-stage LMM modelling, attribution calculation, SATE calculation
?rain_attr

#Bootstrap procedure, and different types of bootstrap
?bootstrap_downwind

#Bootstrap options
?bootstrap_option

#Permutation-based procedure, and different types of permutations
?permutation_ionizer

#Permutation options
?permutation_option

#### Catching user errors ####

# Inconsistency between instrumental prediction variable name vs. its use in downwind_lmm_formula
rain_attr(
  data = oman,
  upwind_lmm_formula = LogRain ~ Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
  #Variable name for storing instrumental prediction from upwind LMM
  instr_pred_name = 'natural_pred',
  instr_pred_type = 'Unconditional',
  #LMM formula fitted to the downwind observations
  downwind_lmm_formula = LogRain ~ natural_pred_typo + Year...2013 + Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation...1km:Target.H.01 + Gauge.Elevation...1km:Target.H.02 + + Gauge.Elevation...1km.1:Target.H.01 + Gauge.Elevation...1km.1:Target.H.02 + (1|TrialDay),
  downwind_logistic_formula = NULL,
  downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
  rain_col_name = 'Rain.Gauge.Measurement',
  upwind_subset = Gauge.Day.Type == 'Upwind' & Year!= 2013,
  downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
  downwind_target_subset = Gauge.Day.Type == 'Target',
  downwind_control_subset = Gauge.Day.Type == 'Control',
  positive_subset = Rain.Gauge.Measurement > 0,
  attr_type = 'ThoEtAl',
  x_downwind_name = c('Year...2013' , 'Year...2014' , 'Year...2016' , 'Year...2017' , 'Year...2018', 'Gauge.Elevation...1km', 'Gauge.Elevation...1km.1', 'natural_pred'),
  target_only = FALSE,
  bootstrap = FALSE,
  bootstrap_option = bootstrap_opt(),
  permutation = FALSE,
  permutation_option = permutation_opt()
)

# Invalid attribution type requested by user
rain_attr(
  data = oman,
  upwind_lmm_formula = LogRain ~ Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
  instr_pred_name = 'natural_pred',
  instr_pred_type = 'Unconditional',
  downwind_lmm_formula = LogRain ~ Year...2013 + Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation...1km:Target.H.01 + Gauge.Elevation...1km:Target.H.02 + + Gauge.Elevation...1km.1:Target.H.01 + Gauge.Elevation...1km.1:Target.H.02 + (1|TrialDay),
  downwind_logistic_formula = NULL,
  downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
  rain_col_name = 'Rain.Gauge.Measurement',
  upwind_subset = Gauge.Day.Type == 'Upwind' & Year!= 2013,
  downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
  downwind_target_subset = Gauge.Day.Type == 'Target',
  downwind_control_subset = Gauge.Day.Type == 'Control',
  positive_subset = Rain.Gauge.Measurement > 0,
  #Invalid attr_type requested
  attr_type = 'Anytype',
  x_downwind_name = c('Year...2013' , 'Year...2014' , 'Year...2016' , 'Year...2017' , 'Year...2018', 'Gauge.Elevation...1km', 'Gauge.Elevation...1km.1', 'natural_pred'),
  target_only = FALSE,
  bootstrap = FALSE,
  bootstrap_option = bootstrap_opt(),
  permutation = FALSE,
  permutation_option = permutation_opt()
)

#### Bootstrap Inference with optional parallelization ####
start_time = Sys.time()
testing_PREB1  = rain_attr(
  data = oman,
  upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
  instr_pred_name = 'natural_pred',
  instr_pred_type = 'Unconditional',
  downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
  downwind_logistic_formula = NULL,
  downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
  rain_col_name = 'Rain.Gauge.Measurement',
  upwind_subset = Gauge.Day.Type == 'Upwind',
  downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
  downwind_target_subset = Gauge.Day.Type == 'Target',
  downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
  attr_type = 'No',
  x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
  target_only = FALSE,
  #TRUE = perform bootstrap inference on attribution and SATE
  bootstrap =T,
  bootstrap_option = bootstrap_opt(
    #Number of bootstrap runs
    B_bootstrap = 500,
    #Type of bootstrap: 'PREB0', 'PREB1', 'PREB2', 'REB0', 'REB1', 'REB2', 'MREB1'
    bootstrap_type = 'PREB1',
    #Whether to bootstrap zeros (i.e., rainfall event indicators)
    bootstrap_zero = F,
    #Probability threshold used for bootstrapping rainfall event indicators
    positive_prob_threshold = NULL,
    #Whether to discretize the bootstrapped rainfall (at raw scale)
    discretize_rain = T,
    #Whether to winsorize individual bootstrapped rainfall to ensure each of them is <= individual_rain_interval[2]
    winsorize_individual_rain = T,
    #lower and upper bounds for adjusting bootstrapped individual rainfall values that are too large when winsorize_individual_rain = TRUE
    individual_rain_interval = c(100, 175),
    #Whether to winsorize the total bootstrapped rainfall (for all downwind observations) to ensure it is within total_rain_interval
    winsorize_total_rain = F,
    # lower and upper bounds for adjusting the total of bootstrapped rainfall values when winsorize_total_rain = TRUE
    total_rain_interval = c(6000, 60000),
    #Seed for reproducibility
    bootstrap_seed = 123,
    #Whether to parallelize the bootstraps
    bootstrap_parallel = T,
    #Number of workers for parallelization
    bootstrap_parallel_num_worker = parallel::detectCores() - 1,
    #Confidence level of bootstrap CI
    CI_level = 0.95
  ),
  permutation = FALSE,
  permutation_option = permutation_opt()
)
end_time = Sys.time()
end_time - start_time

#Bootstrap CI result for attribution and SATE, as well as the parameters of all fitted models e.g., downwind lmm model, downwind propensity score model
testing_PREB1$bootstrap_CI_result

#Bootstrap p-value result (proportion of bootstrapped estimates that are < 0) for attribution and SATE, as well as the parameters of all fitted models (but less meaningful for these parameters)
testing_PREB1$bootstrap_p_value_result

#Bootstrap distributions for attribution and SATE, as well as the parameters of all fitted models, and the bootstrapped responses (could include offset term) of the downwind_lmm_formula
head(testing_PREB1$bootstrap_result$hatattr)
head(testing_PREB1$bootstrap_result$hatsate)

#Plots of bootstrap distribution for attribution and SATE
ggpubr::ggarrange(plotlist =testing_PREB1$bootstrap_plot_result$hatattr)
ggpubr::ggarrange(plotlist = testing_PREB1$bootstrap_plot_result$hatsate)



#### Permutation Inference ####
start_time = Sys.time()
testing_perm = rain_attr(
  data = oman,
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
  attr_type = 'No',
  x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
  target_only = FALSE,
  bootstrap =F,
  bootstrap_option = NULL,
  permutation = T,
  permutation_option = permutation_opt(
    # Number of permutation run
    B_permutation = 500,
    # Whether to permute the elements within each row of ionizer_operation_input
    permute_between_ionizer = T,
    # Whether to permute the rows of ionizer_operation_input
    permute_all_ionizers_between_day = F,
    # Whether to permute the rows of gauge-day level of ionizer_operation_input
    permute_between_gaugeday = T,
    # Data frame containing day-level information on the operating states (1 for active, 0 for inactive) of all ionizers for all days of the trial
    ionizer_operation_input = ionizer_operation,
    # A matrix containing gauge-day (row) level information on whether each gauge is downwind of an ionizer (column) on a particular day
    gaugeday_downwind_input = gaugeday_downwind,
    # Vector of column names in the supplied data, that corresponds to each ionizer's target indicator
    data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
    # Column name in ionizer_operation_input capturing year information
    ionizer_operation_year_column_name = 'Year',
    # Column name in ionizer_operation_input capturing day information
    ionizer_operation_day_column_name = 'TrialDay',
    #Seed for reproducibility
    permutation_seed = 123,
    #Whether to parallelize the permutations
    permutation_parallel = T,
    #Number of workers for parallelization
    permutation_parallel_num_worker = parallel::detectCores() - 1
  )
)
end_time = Sys.time()
end_time - start_time

#### Bootstrap and Permutation Inference ####
start_time = Sys.time()
test  = rain_attr(
  data = oman,
  upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
  instr_pred_name = 'natural_pred',
  instr_pred_type = 'Unconditional',
  downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
  downwind_logistic_formula = NULL,
  downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
  rain_col_name = 'Rain.Gauge.Measurement',
  upwind_subset = Gauge.Day.Type == 'Upwind',
  downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
  downwind_target_subset = Gauge.Day.Type == 'Target',
  downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
  attr_type = 'No',
  x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
  target_only = FALSE,
  #TRUE = perform bootstrap inference on attribution and SATE
  bootstrap =T,
  bootstrap_option = bootstrap_opt(
    #Number of bootstrap runs
    B_bootstrap = 500,
    #Type of bootstrap: 'PREB0', 'PREB1', 'PREB2', 'REB0', 'REB1', 'REB2', 'MREB1'
    bootstrap_type = 'PREB1',
    #Whether to bootstrap zeros (i.e., rainfall event indicators)
    bootstrap_zero = F,
    #Probability threshold used for bootstrapping rainfall event indicators
    positive_prob_threshold = NULL,
    #Whether to discretize the bootstrapped rainfall (at raw scale)
    discretize_rain = T,
    #Whether to winsorize individual bootstrapped rainfall to ensure each of them is <= individual_rain_interval[2]
    winsorize_individual_rain = T,
    #lower and upper bounds for adjusting bootstrapped individual rainfall values that are too large when winsorize_individual_rain = TRUE
    individual_rain_interval = c(100, 175),
    #Whether to winsorize the total bootstrapped rainfall (for all downwind observations) to ensure it is within total_rain_interval
    winsorize_total_rain = F,
    # lower and upper bounds for adjusting the total of bootstrapped rainfall values when winsorize_total_rain = TRUE
    total_rain_interval = c(6000, 60000),
    #Seed for reproducibility
    bootstrap_seed = 123,
    #Whether to parallelize the bootstraps
    bootstrap_parallel = T,
    #Number of workers for parallelization
    bootstrap_parallel_num_worker = parallel::detectCores() - 1,
    #Confidence level of bootstrap CI
    CI_level = 0.95
  ),
  permutation = T,
  permutation_option = permutation_opt(
    # Number of permutation run
    B_permutation = 500,
    # Whether to permute the elements within each row of ionizer_operation_input
    permute_between_ionizer = T,
    # Whether to permute the rows of ionizer_operation_input
    permute_all_ionizers_between_day = F,
    # Whether to permute the rows of gauge-day level of ionizer_operation_input
    permute_between_gaugeday = T,
    # Data frame containing day-level information on the operating states (1 for active, 0 for inactive) of all ionizers for all days of the trial
    ionizer_operation_input = ionizer_operation,
    # A matrix containing gauge-day (row) level information on whether each gauge is downwind of an ionizer (column) on a particular day
    gaugeday_downwind_input = gaugeday_downwind,
    # Vector of column names in the supplied data, that corresponds to each ionizer's target indicator
    data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
    # Column name in ionizer_operation_input capturing year information
    ionizer_operation_year_column_name = 'Year',
    # Column name in ionizer_operation_input capturing day information
    ionizer_operation_day_column_name = 'TrialDay',
    #Seed for reproducibility
    permutation_seed = 123,
    #Whether to parallelize the permutations
    permutation_parallel = T,
    #Number of workers for parallelization
    permutation_parallel_num_worker = parallel::detectCores() - 1
  )
)
end_time = Sys.time()
end_time - start_time


#Permutation p-value (proportion of permuted estimates that are greater than or equal to original estimate) for attribution and SATE
testing_perm$permutation_p_value_result

#Permutation distribution for attribution and SATE
testing_perm$permutation_result

#Plots of permutation distribution for attribution and SATE
ggpubr::ggarrange(plotlist = testing_perm$permutation_plot_result$hatattr)
ggpubr::ggarrange(plotlist = testing_perm$permutation_plot_result$hatsate)


#### EDA visualizations ####

#"num_obs_days", "num_obs_days_by_year", "hist_day_group_sizes", "qq_rain", "ts_by_type", "ts_by_gauge", "ts_by_gauge_interactive", "map_static", "map_dynamic"

test_eda = eda(eda_type = "ts_by_gauge_interactive",
              data = oman,
              rain_col_name = 'Rain.Gauge.Measurement',
              day_column_name = 'TrialDay',
              year_column_name = 'Year',
              use_raw = F,
              gauge_id_column_name = 'Gauge.ID',
              ts_focus_gauge = c(1,2),
              longlat_column_names = c("Gauge.Longitude", "Gauge.Latitude"),
              long_lim = c(55,60),
              lat_lim = c(22,25),
              input_sf = rnaturalearth::ne_countries(scale = "large", country = "Oman", returnclass = "sf"),
              # input_sf = NULL,
              ionizer_location_df = ionizer_location,
              ionizer_id_column_name = 'Ionizer',
              ionizer_longlat_column_names = c("Longitude","Latitude"),
              elev_contour = TRUE,
              elev_resolution = 2,
              focus_year = c(2013,2014),
              fps = 10,
              upwind_subset = Gauge.Day.Type == 'Upwind',
              downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
              downwind_target_subset = Gauge.Day.Type == 'Target',
              downwind_control_subset = Gauge.Day.Type == 'Control',
              positive_subset = Rain.Gauge.Measurement > 0
)
test_eda

high_gauges = unique(oman$Gauge.ID[oman$Gauge.Elevation > median(oman$Gauge.Elevation)])

test_focus =eda(eda_type = "ts_by_gauge",
                data = oman,
                rain_col_name = 'Rain.Gauge.Measurement',
                day_column_name = 'TrialDay',
                year_column_name = 'Year',
                use_raw = F,
                gauge_id_column_name = 'Gauge.ID',
                ts_focus_gauge = high_gauges,
                longlat_column_names = c("Gauge.Longitude", "Gauge.Latitude"),
                long_lim = c(55,60),
                lat_lim = c(22,25),
                input_sf = rnaturalearth::ne_countries(scale = "large", country = "Oman", returnclass = "sf"),
                # input_sf = NULL,
                ionizer_location_df = ionizer_location,
                ionizer_id_column_name = 'Ionizer',
                ionizer_longlat_column_names = c("Longitude","Latitude"),
                focus_year = c(2013,2014),
                fps = 10,
                upwind_subset = Gauge.Day.Type == 'Upwind',
                downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
                downwind_target_subset = Gauge.Day.Type == 'Target',
                downwind_control_subset = Gauge.Day.Type == 'Control',
                positive_subset = Rain.Gauge.Measurement > 0
) + ggplot2::guides(color = 'none')

test_focus


animation_out = eda(eda_type = "map_dynamic",
                    data = oman,
                    rain_col_name = 'Rain.Gauge.Measurement',
                    day_column_name = 'TrialDay',
                    year_column_name = 'Year',
                    use_raw = F,
                    gauge_id_column_name = 'Gauge.ID',
                    ts_focus_gauge = c(1,2),
                    longlat_column_names = c("Gauge.Longitude", "Gauge.Latitude"),
                    long_lim = c(55,60),
                    lat_lim = c(22,25),
                    input_sf = rnaturalearth::ne_countries(scale = "large", country = "Oman", returnclass = "sf"),
                    # input_sf = NULL,
                    ionizer_location_df = ionizer_location,
                    # ionizer_location_df = NULL,
                    ionizer_id_column_name = 'Ionizer',
                    ionizer_longlat_column_names = c("Longitude","Latitude"),
                    elev_contour = FALSE,
                    elev_resolution = 2,
                    focus_year = c(2013),
                    fps = 30,
                    upwind_subset = Gauge.Day.Type == 'Upwind',
                    downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
                    downwind_target_subset = Gauge.Day.Type == 'Target',
                    downwind_control_subset = Gauge.Day.Type == 'Control',
                    positive_subset = Rain.Gauge.Measurement > 0
)

gganimate::animate(
  animation_out,
  nframes = 740,            # total frames
  fps = 30,
  width = 800,
  height = 600,
  res = 100,
  renderer = gganimate::file_renderer(
    dir = "frames",
    prefix = "frame",
    overwrite = TRUE
  )
)


#For checking: apo should be 0.111206, apl should be 0.1251201 for attr_type = 'ChambersEtAl'
# apo  = 0.06265952, apl =  0.0668482 for attr_type = 'ThoEtAl'
# # > asd$hatsate$sate.mb; asd$hatsate$sate.ipw; asd$hatsate$sate.ipw.l
# [1] 0.1143799
# [1] 0.07401868
# [1] 0.1122864
#
# > asd$hatsate$sate.ipw.ma; asd$hatsate$sate.aipw
# [1] 0.06543017
# [1] 0.07736621

#
#Upwind:
#> asd$fitted_models$upwind_lmm_fit
# Linear mixed model fit by REML ['lmerMod']
# Formula: LogRain ~ Gauge.Elevation + Steering.Wind.Speed + Total.Totals +      PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure +      (1 | TrialDay)
# Data: data[upwind & positive, ]
# REML criterion at convergence: 5345.72
# Random effects:
#   Groups   Name        Std.Dev.
# TrialDay (Intercept) 0.6449
# Residual             1.2650
# Number of obs: 1545, groups:  TrialDay, 292
# Fixed Effects:
#   (Intercept)            Gauge.Elevation        Steering.Wind.Speed               Total.Totals        PC2.Dry.Temperature      PC1.Relative.Humidity  PC1.Ground.Level.Pressure
# -1.43966                    0.43401                   -0.09641                    0.03268                    0.14478                    0.17787                   -0.05232


#> asd$fitted_models$downwind_lmm_fit
# Linear mixed model fit by REML ['lmerMod']
# Formula: LogRain ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 +      Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 +
#   Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 +      Gauge.Elevation:Target.H.02 + (1 | TrialDay)
# Data: data[downwind & positive, ]
# REML criterion at convergence: 14764.41
# Random effects:
#   Groups   Name        Std.Dev.
# TrialDay (Intercept) 0.5232
# Residual             1.3623
# Number of obs: 4168, groups:  TrialDay, 488
# Fixed Effects:
#   (Intercept)              Gauge.Elevation                 natural_pred                  Target.H.01                  Target.H.02                  Target.H.03
# 0.27975                     -0.12503                      0.85569                      0.28922                      0.25764                      0.23818
# Target.H.04                  Target.H.05                  Target.H.06                  Target.H.07                  Target.H.08                  Target.H.09
# -0.15266                      0.43253                     -0.20325                      0.22201                      0.05850                      0.48147
# Target.H.10  Gauge.Elevation:Target.H.01  Gauge.Elevation:Target.H.02
# 0.03444                     -0.08662                     -0.21665


# data = oman
# upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay)
# instr_pred_name = 'natural_pred'
# downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay)
# downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02
# downwind_propensity_formula = Gauge.Day.Type == 'Target' ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure
# rain_col_name = 'Rain.Gauge.Measurement'
# #
#
# asd  = rain_attr(data = oman,
#                  upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                  instr_pred_name = 'natural_pred',
#                  instr_pred_type = 'Unconditional',
#                  downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                  downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                  downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                  rain_col_name = 'Rain.Gauge.Measurement',
#                  upwind_subset = Gauge.Day.Type == 'Upwind',
#                  downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                  downwind_target_subset = Gauge.Day.Type == 'Target',
#                  downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                  attr_type = 'ThoEtAl',
#                  x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                  target_only = FALSE)
#
# asd$hatsate$sate.mb; asd$hatsate$sate.ipw; asd$hatsate$sate.ipw.l
# asd$hatsate$sate.ipw.ma; asd$hatsate$sate.aipw
# asd$hatattr$apl; asd$hatattr$apo
#
#
# #To replicate Table 6 of JRSSA
# oman_in = oman
# replicate_table6 = rain_attr(data = oman_in,
#                              upwind_lmm_formula = LogRain ~ Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                              instr_pred_name = 'natural_pred',
#                              instr_pred_type = 'Unconditional',
#                              downwind_lmm_formula = LogRain ~ Year...2013 + Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation...1km:Target.H.01 + Gauge.Elevation...1km:Target.H.02 + + Gauge.Elevation...1km.1:Target.H.01 + Gauge.Elevation...1km.1:Target.H.02 + (1|TrialDay),
#                              downwind_logistic_formula = NULL,
#                              downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                              rain_col_name = 'Rain.Gauge.Measurement',
#                              upwind_subset = Gauge.Day.Type == 'Upwind' & Year!= 2013,
#                              downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                              downwind_target_subset = Gauge.Day.Type == 'Target',
#                              downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                              attr_type = 'ThoEtAl',
#                              x_downwind_name = c('Year...2013' , 'Year...2014' , 'Year...2016' , 'Year...2017' , 'Year...2018', 'Gauge.Elevation...1km', 'Gauge.Elevation...1km.1', 'natural_pred'),
#                              target_only = FALSE)
#
#
# #Table 1 of JRSSA
# replicate_table6$all_fitted_models$downwind_propensity_fit
#
# #Table 2 of JRSSA
# lme4::fixef(replicate_table6$all_fitted_models$downwind_lmm_fit)
#
# #Table 3 of JRSSA
# as.data.frame(lme4::VarCorr(replicate_table6$all_fitted_models$downwind_lmm_fit))
#
# #Table 4 of JRSSA
# replicate_table6$all_fitted_models$downwind_positive_target_lmm_fit
# replicate_table6$all_fitted_models$downwind_positive_control_lmm_fit
#
# #Table 5
# as.data.frame(lme4::VarCorr(replicate_table6$all_fitted_models$downwind_positive_target_lmm_fit))
# as.data.frame(lme4::VarCorr(replicate_table6$all_fitted_models$downwind_positive_control_lmm_fit))
#
# #DONE: Check why our AIPW estimate is different from Table 6 in JRSSA which gives 0.073 but here we get 0.07583943
# #ANS: The reason is probably because Ray was computing sate.aipw using the incorrect expression (relationship between AIPW and MB) below equation (7) instead of directly using his equation (7).
# #     This can be verified since our wrong.sate.aipw (computed using the incorrect relationship between AIPW and MB below equation (7) ) gives 0.07288667 which is same as the reported 0.073 of AIPW in Table 6
#
# #Estimate row for Table 6 of LogRain
# unlist(replicate_table6$hatsate)
#
#
# lme4::fixef(replicate_table6$all_fitted_models$upwind_lmm_fit)

# B_bootstrap = 3
# bootstrap_type = 'REB0'
# bootstrap_zero = TRUE
# positive_prob_threshold = NULL
# discretize_rain = TRUE
# winsorize_individual_rain = TRUE
# winsorize_total_rain = TRUE
#
# ori_data = asd$data
#
# #upwind = oman$Gauge.Day.Type == 'Upwind'
# downwind = oman$Gauge.Day.Type  %in% c('Target','Control')
# ori_positive = oman$Rain.Gauge.Measurement > 0
# rain_col_name = 'Rain.Gauge.Measurement'
# x_downwind_name = c('Gauge.Elevation', 'natural_pred')
# ori_attr_est = asd$hatattr
# ori_sate_est = asd$hatsate
# ori_fitted_models = asd$fitted_models
#
# #
# z_downwind_name = setdiff(names(lme4::fixef(asd$fitted_models$downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
# downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)
# target_only = F
# attr_type = 'ThoEtAl'

# asd2  = rain_attr(data = oman,
#                   upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                   instr_pred_name = 'natural_pred',
#                   instr_pred_type = 'Unconditional',
#                   downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                   downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                   downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                   rain_col_name = 'Rain.Gauge.Measurement',
#                   upwind_subset = Gauge.Day.Type == 'Upwind',
#                   downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                   downwind_target_subset = Gauge.Day.Type == 'Target',
#                   downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                   attr_type = 'ThoEtAl',
#                   x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                   target_only = FALSE,
#                   bootstrap =T,
#                   bootstrap_option = bootstrap_opt(B_bootstrap = 3,
#                                                       bootstrap_type = 'REB1',
#                                                       bootstrap_zero = T,
#                                                       positive_prob_threshold = NULL,
#                                                       discretize_rain = T,
#                                                       winsorize_individual_rain = T,
#                                                       winsorize_total_rain = T,
#                                                       CI_level = 0.95
#                   )
# )
# c(asd2$hatattr$apo, asd2$hatattr$apl)
# c(asd$hatattr$apo, asd$hatattr$apl)
# unlist(asd2$hatsate)
# unlist(asd$hatsate)
#
# asd2$bootstrap_result$hatattr
# asd2$bootstrap_result$hatsate
# asd2$bootstrap_result$downwind_lmm_param
# asd2$bootstrap_result$downwind_logistic_param
# asd2$bootstrap_result$downwind_propensity_param
# asd2$bootstrap_result$downwind_positive_target_lmm_param
# asd2$bootstrap_result$downwind_positive_control_lmm_param
# asd2$bootstrap_result$downwind_LogRain[1:3,1:10]
# apply(asd2$bootstrap_result$downwind_LogRain, 1 , function(x){sum(!is.na(x))})
#
#Trying to replicate PREB-1 bootstrap paper results:
# set.seed(123)
# asd3_PREB1  = rain_attr(data = oman,
#                         upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                         instr_pred_name = 'natural_pred',
#                         instr_pred_type = 'Unconditional',
#                         downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                         downwind_logistic_formula = NULL,
#                         downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                         rain_col_name = 'Rain.Gauge.Measurement',
#                         upwind_subset = Gauge.Day.Type == 'Upwind',
#                         downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                         downwind_target_subset = Gauge.Day.Type == 'Target',
#                         downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                         attr_type = 'No',
#                         x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                         target_only = FALSE,
#                         bootstrap =T,
#                         bootstrap_option = bootstrap_opt(B_bootstrap = 3,
#                                                             bootstrap_type = 'PREB1',
#                                                             bootstrap_zero = F,
#                                                             positive_prob_threshold = NULL,
#                                                             discretize_rain = F,
#                                                             winsorize_individual_rain = F,
#                                                             winsorize_total_rain = F,
#                                                             CI_level = 0.95
#                         )
# )
# #Verified to be equivalent to the first 3 bootstrap runs of bootstrap paper for PREB1
# #Only has minor difference in terms of 4th or 5th decimal points for some results, which could be due to the use of lme4::lmer() but we were using nlme::lme() in the bootstrap paper
# c(asd3_PREB1$hatattr$apo,asd3_PREB1$hatattr$apl, asd3_PREB1$hatsate$sate.mb, asd3_PREB1$hatsate$sate.ipw, asd3_PREB1$hatsate$sate.ipw.l)
# t(cbind(asd3_PREB1$bootstrap_result$hatattr,asd3_PREB1$bootstrap_result$hatsate))
# asd3_PREB1$bootstrap_result$downwind_lmm_param
#
# load("D:/Postdoc/Bootstrap Paper/R Codes/Rdata/real_data_B500.rda")
# PREB1.quantities.bootstrap.distribution[,1:3] - t(cbind(asd3_PREB1$bootstrap_result$hatattr*100,asd3_PREB1$bootstrap_result$hatsate[,1:3]))
# max(abs(PREB1.quantities.bootstrap.distribution[,1:3] - t(cbind(asd3_PREB1$bootstrap_result$hatattr*100,asd3_PREB1$bootstrap_result$hatsate[,1:3]))))
#
# PREB1.result$PREB1.result[1:3,-18] - asd3_PREB1$bootstrap_result$downwind_lmm_param
# max(abs(PREB1.result$PREB1.result[1:3,-18] - asd3_PREB1$bootstrap_result$downwind_lmm_param))

#Testing REB2
# set.seed(123)
# asd3_PREB2  = rain_attr(data = oman,
#                         upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                         instr_pred_name = 'natural_pred',
#                         instr_pred_type = 'Unconditional',
#                         downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                         downwind_logistic_formula = NULL,
#                         downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                         rain_col_name = 'Rain.Gauge.Measurement',
#                         upwind_subset = Gauge.Day.Type == 'Upwind',
#                         downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                         downwind_target_subset = Gauge.Day.Type == 'Target',
#                         downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                         attr_type = 'No',
#                         x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                         target_only = FALSE,
#                         bootstrap =T,
#                         bootstrap_option = bootstrap_opt(B_bootstrap = 10,
#                                                             bootstrap_type = 'PREB2',
#                                                             bootstrap_zero = F,
#                                                             positive_prob_threshold = NULL,
#                                                             discretize_rain = F,
#                                                             winsorize_individual_rain = F,
#                                                             winsorize_total_rain = F,
#                                                             CI_level = 0.95
#                         )
# )
# apply(asd3_PREB2$bootstrap_result$downwind_positive_target_lmm_param,2,mean)
# asd3_PREB2$all_fitted_models$downwind_positive_target_lmm_fit
#
#
# #Testing REB2
# set.seed(123)
# asd3_REB2  = rain_attr(data = oman,
#                        upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                        instr_pred_name = 'natural_pred',
#                        instr_pred_type = 'Unconditional',
#                        downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                        downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                        downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                        rain_col_name = 'Rain.Gauge.Measurement',
#                        upwind_subset = Gauge.Day.Type == 'Upwind',
#                        downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                        downwind_target_subset = Gauge.Day.Type == 'Target',
#                        downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                        attr_type = 'No',
#                        x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                        target_only = FALSE,
#                        bootstrap =T,
#                        bootstrap_option = bootstrap_opt(B_bootstrap = 10,
#                                                            bootstrap_type = 'REB2',
#                                                            bootstrap_zero = F,
#                                                            positive_prob_threshold = NULL,
#                                                            discretize_rain = F,
#                                                            winsorize_individual_rain = F,
#                                                            winsorize_total_rain = F,
#                                                            CI_level = 0.95
#                        )
# )
# apply(asd3_REB2$bootstrap_result$downwind_positive_target_lmm_param,2,mean)
# asd3_REB2$all_fitted_models$downwind_positive_target_lmm_fit
#
# #Testing bootstrap_zero
# set.seed(123)
# asd3_REB2_bootstrap_zero  = rain_attr(data = oman,
#                                       upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                                       instr_pred_name = 'natural_pred',
#                                       instr_pred_type = 'Unconditional',
#                                       downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                                       downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                                       downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                                       rain_col_name = 'Rain.Gauge.Measurement',
#                                       upwind_subset = Gauge.Day.Type == 'Upwind',
#                                       downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                                       downwind_target_subset = Gauge.Day.Type == 'Target',
#                                       downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                                       attr_type = 'No',
#                                       x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                       target_only = FALSE,
#                                       bootstrap =T,
#                                       bootstrap_option = bootstrap_opt(B_bootstrap = 10,
#                                                                           bootstrap_type = 'REB2',
#                                                                           bootstrap_zero = T,
#                                                                           positive_prob_threshold = NULL,
#                                                                           discretize_rain = F,
#                                                                           winsorize_individual_rain = F,
#                                                                           winsorize_total_rain = F,
#                                                                           CI_level = 0.95
#                                       )
# )
# apply(asd3_REB2_bootstrap_zero$bootstrap_result$downwind_positive_target_lmm_param,2,mean)
# asd3_REB2_bootstrap_zero$all_fitted_models$downwind_positive_target_lmm_fit
# asd3_REB2_bootstrap_zero$bootstrap_CI_result
# asd3_REB2_bootstrap_zero$bootstrap_p_value_result
# ggpubr::ggarrange(plotlist =asd3_REB2_bootstrap_zero$bootstrap_plot_result$hatattr)
# ggpubr::ggarrange(plotlist = asd3_REB2_bootstrap_zero$bootstrap_plot_result$hatsate)
#
# #Testing bootstrap_zero and discretize
# set.seed(123)
# asd3_REB2_bootstrap_zero_discrete  = rain_attr(data = oman,
#                                                upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                                                instr_pred_name = 'natural_pred',
#                                                instr_pred_type = 'Unconditional',
#                                                downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                                                downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                                                downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                                                rain_col_name = 'Rain.Gauge.Measurement',
#                                                upwind_subset = Gauge.Day.Type == 'Upwind',
#                                                downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                                                downwind_target_subset = Gauge.Day.Type == 'Target',
#                                                downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                                                attr_type = 'No',
#                                                x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                                target_only = FALSE,
#                                                bootstrap =T,
#                                                bootstrap_option = bootstrap_opt(B_bootstrap = 10,
#                                                                                    bootstrap_type = 'REB2',
#                                                                                    bootstrap_zero = T,
#                                                                                    positive_prob_threshold = NULL,
#                                                                                    discretize_rain = T,
#                                                                                    winsorize_individual_rain = T,
#                                                                                    winsorize_total_rain = T,
#                                                                                    CI_level = 0.95
#                                                )
# )
# apply(asd3_REB2_bootstrap_zero_discrete$bootstrap_result$downwind_positive_target_lmm_param,2,mean)
# asd3_REB2_bootstrap_zero_discrete$all_fitted_models$downwind_positive_target_lmm_fit
# asd3_REB2_bootstrap_zero_discrete$bootstrap_CI_result
# asd3_REB2_bootstrap_zero_discrete$bootstrap_p_value_result
# ggpubr::ggarrange(plotlist =asd3_REB2_bootstrap_zero_discrete$bootstrap_plot_result$hatattr)
# ggpubr::ggarrange(plotlist = asd3_REB2_bootstrap_zero_discrete$bootstrap_plot_result$hatsate)
#
# asd3_REB2_bootstrap_zero_discrete$bootstrap_result$downwind_positive_target_lmm_param
# asd3_REB2_bootstrap_zero$bootstrap_result$downwind_positive_target_lmm_param
# asd3_REB2$bootstrap_result$downwind_positive_target_lmm_param
# asd3_PREB2$bootstrap_result$downwind_positive_target_lmm_param

#TODO: verify the results with bootstrap paper for REB1 and MREB1 (need to modify the prob argument in sample() to match bootstrap paper), and then with
#D:\Postdoc\Simulation\Replicate ISR Results\Bootstrap Analysis with generate_zero_T and scaled_h_sampling and Correct Scaling REB1 using Oman Data.R
#particularly for bootstrap_zero = T, discretize_rain = T, winsorize_individual_rain = T, winsorize_total_rain = T - but this might need some restructuring of the bootstrap_downwind() since the ordering of sampling matters!

#Testing permutation
# set.seed(123)
# testing = rain_attr(data = oman,
#                     upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                     instr_pred_name = 'natural_pred',
#                     instr_pred_type = 'Unconditional',
#                     downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                     downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                     downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                     rain_col_name = 'Rain.Gauge.Measurement',
#                     upwind_subset = Gauge.Day.Type == 'Upwind',
#                     downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                     downwind_target_subset = Gauge.Day.Type == 'Target',
#                     downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                     attr_type = 'No',
#                     x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                     target_only = FALSE,
#                     bootstrap =F,
#                     bootstrap_option = NULL,
#                     permutation = T,
#                     permutation_option = permutation_opt(
#                       B_permutation = 5,
#                       permute_between_ionizer = T,
#                       permute_all_ionizers_between_day = T,
#                       permute_between_gaugeday = T,
#                       ionizer_operation_input = ionizer_operation,
#                       gaugeday_downwind_input = gaugeday_downwind,
#                       year_ionizer_list =
#                         list(
#                           '2013' = c('H1','H2'),
#                           '2014' = c('H1','H2','H3','H4'),
#                           '2015' = c('H1','H2','H3','H4','H5','H6'),
#                           '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
#                           '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
#                           '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
#                         ),
#                       data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
#                       ionizer_operation_year_column_name = 'Year',
#                       ionizer_operation_day_column_name = 'TrialDay'
#                     )
# )
# testing$permutation_result
# testing$permutation_p_value_result
# ggpubr::ggarrange(plotlist = testing$permutation_plot_result$hatattr)
# ggpubr::ggarrange(plotlist = testing$permutation_plot_result$hatsate)


#Replicating previous permutation analysis for permute_between_gaugeday = F
#Note we need to use sample.kind = 'Rounding' due to previous analysis loaded RData8.Rdata, which caused sample.kind = 'Rounding' from older R version instead of sample.kind = 'Rejection' in the latest R version
# RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
# set.seed(123)
# my_perm_result_TT_ChambersEtAl = rain_attr(data = oman,
#                                            upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                                            instr_pred_name = 'natural_pred',
#                                            instr_pred_type = 'Unconditional',
#                                            downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                                            downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                                            downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                                            rain_col_name = 'Rain.Gauge.Measurement',
#                                            upwind_subset = Gauge.Day.Type == 'Upwind',
#                                            downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                                            downwind_target_subset = Gauge.Day.Type == 'Target',
#                                            downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                                            attr_type = 'ChambersEtAl',
#                                            x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                            target_only = FALSE,
#                                            bootstrap =F,
#                                            bootstrap_option = NULL,
#                                            permutation = T,
#                                            permutation_option = permutation_opt(
#                                              B_permutation = 6,
#                                              permute_between_ionizer = T,
#                                              permute_all_ionizers_between_day = T,
#                                              permute_between_gaugeday = F,
#                                              ionizer_operation_input = ionizer_operation,
#                                              gaugeday_downwind_input = gaugeday_downwind,
#                                              year_ionizer_list =
#                                                list(
#                                                  '2013' = c('H1','H2'),
#                                                  '2014' = c('H1','H2','H3','H4'),
#                                                  '2015' = c('H1','H2','H3','H4','H5','H6'),
#                                                  '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
#                                                  '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
#                                                  '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
#                                                ),
#                                              data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
#                                              ionizer_operation_year_column_name = 'Year',
#                                              ionizer_operation_day_column_name = 'TrialDay'
#                                            )
# )
#
# load('D:/Postdoc/Simulation/Replicate ISR Results/Rdata/permutation_result_Oman_Trial_Data_perm_row_between_gauge_day_F.Rdata')
# max(abs(perm_result_TT$perm_attribution_ChambersEtAl_matrix[1:6,c('apo','apl')] - my_perm_result_TT_ChambersEtAl$permutation_result$hatattr))
#
#
# RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
# set.seed(123)
# my_perm_result_TT_ThoEtAl = rain_attr(data = oman,
#                                            upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                                            instr_pred_name = 'natural_pred',
#                                            instr_pred_type = 'Unconditional',
#                                            downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                                            downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                                            downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                                            rain_col_name = 'Rain.Gauge.Measurement',
#                                            upwind_subset = Gauge.Day.Type == 'Upwind',
#                                            downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                                            downwind_target_subset = Gauge.Day.Type == 'Target',
#                                            downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                                            attr_type = 'ThoEtAl',
#                                            x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                            target_only = FALSE,
#                                            bootstrap =F,
#                                            bootstrap_option = NULL,
#                                            permutation = T,
#                                            permutation_option = permutation_opt(
#                                              B_permutation = 6,
#                                              permute_between_ionizer = T,
#                                              permute_all_ionizers_between_day = T,
#                                              permute_between_gaugeday = F,
#                                              ionizer_operation_input = ionizer_operation,
#                                              gaugeday_downwind_input = gaugeday_downwind,
#                                              year_ionizer_list =
#                                                list(
#                                                  '2013' = c('H1','H2'),
#                                                  '2014' = c('H1','H2','H3','H4'),
#                                                  '2015' = c('H1','H2','H3','H4','H5','H6'),
#                                                  '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
#                                                  '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
#                                                  '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
#                                                ),
#                                              data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
#                                              ionizer_operation_year_column_name = 'Year',
#                                              ionizer_operation_day_column_name = 'TrialDay'
#                                            )
# )
# max(abs(perm_result_TT$perm_attribution_proposed_matrix[1:6,c('apo','apl')] - my_perm_result_TT_ThoEtAl$permutation_result$hatattr))
#
#
# #Replicate previosu analysis with permute_between_gaugeday = T
# RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
# set.seed(123)
# my_perm_result_TT_ChambersEtAl = rain_attr(data = oman,
#                                            upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                                            instr_pred_name = 'natural_pred',
#                                            instr_pred_type = 'Unconditional',
#                                            downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                                            downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                                            downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                                            rain_col_name = 'Rain.Gauge.Measurement',
#                                            upwind_subset = Gauge.Day.Type == 'Upwind',
#                                            downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                                            downwind_target_subset = Gauge.Day.Type == 'Target',
#                                            downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                                            attr_type = 'ChambersEtAl',
#                                            x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                            target_only = FALSE,
#                                            bootstrap =F,
#                                            bootstrap_option = NULL,
#                                            permutation = T,
#                                            permutation_option = permutation_opt(
#                                              B_permutation = 6,
#                                              permute_between_ionizer = T,
#                                              permute_all_ionizers_between_day = T,
#                                              permute_between_gaugeday = T,
#                                              ionizer_operation_input = ionizer_operation,
#                                              gaugeday_downwind_input = gaugeday_downwind,
#                                              year_ionizer_list =
#                                                list(
#                                                  '2013' = c('H1','H2'),
#                                                  '2014' = c('H1','H2','H3','H4'),
#                                                  '2015' = c('H1','H2','H3','H4','H5','H6'),
#                                                  '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
#                                                  '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
#                                                  '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
#                                                ),
#                                              data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
#                                              ionizer_operation_year_column_name = 'Year',
#                                              ionizer_operation_day_column_name = 'TrialDay'
#                                            )
# )
#
# load('D:/Postdoc/Simulation/Replicate ISR Results/Rdata/permutation_result_Oman_Trial_Data_perm_row_between_gauge_day_T.Rdata')
# max(abs(perm_result_TT$perm_attribution_ChambersEtAl_matrix[1:6,c('apo','apl')] - my_perm_result_TT_ChambersEtAl$permutation_result$hatattr))
#
# RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
# set.seed(123)
# my_perm_result_TT_ThoEtAl = rain_attr(data = oman,
#                                        upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                                        instr_pred_name = 'natural_pred',
#                                        instr_pred_type = 'Unconditional',
#                                        downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                                        downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                                        downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                                        rain_col_name = 'Rain.Gauge.Measurement',
#                                        upwind_subset = Gauge.Day.Type == 'Upwind',
#                                        downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                                        downwind_target_subset = Gauge.Day.Type == 'Target',
#                                        downwind_control_subset = Gauge.Day.Type == 'Control', positive_subset = Rain.Gauge.Measurement > 0,
#                                        attr_type = 'ThoEtAl',
#                                        x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                        target_only = FALSE,
#                                        bootstrap =F,
#                                        bootstrap_option = NULL,
#                                        permutation = T,
#                                        permutation_option = permutation_opt(
#                                          B_permutation = 6,
#                                          permute_between_ionizer = T,
#                                          permute_all_ionizers_between_day = T,
#                                          permute_between_gaugeday = T,
#                                          ionizer_operation_input = ionizer_operation,
#                                          gaugeday_downwind_input = gaugeday_downwind,
#                                          year_ionizer_list =
#                                            list(
#                                              '2013' = c('H1','H2'),
#                                              '2014' = c('H1','H2','H3','H4'),
#                                              '2015' = c('H1','H2','H3','H4','H5','H6'),
#                                              '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
#                                              '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
#                                              '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
#                                            ),
#                                          data_target_column_names = c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10"),
#                                          ionizer_operation_year_column_name = 'Year',
#                                          ionizer_operation_day_column_name = 'TrialDay'
#                                        )
# )
# max(abs(perm_result_TT$perm_attribution_proposed_matrix[1:6,c('apo','apl')] - my_perm_result_TT_ThoEtAl$permutation_result$hatattr))

# RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

#
# set.seed(123)
# sample(1:2)
#
# set.seed(123)
# sample(c('a','b'))

#TODO: Do a full-scale replication of previous permutation analysis, bootstrap analysis in Bootstrap paper, and previous bootstrap analysis (involving bootstrap_zero) to verify the correctness of all functions

