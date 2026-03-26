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
    # A list supplying information on deployed ionizers for each year
    year_ionizer_list =
      list(
        '2013' = c('H1','H2'),
        '2014' = c('H1','H2','H3','H4'),
        '2015' = c('H1','H2','H3','H4','H5','H6'),
        '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
        '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
        '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
      ),
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
    # A list supplying information on deployed ionizers for each year
    year_ionizer_list =
      list(
        '2013' = c('H1','H2'),
        '2014' = c('H1','H2','H3','H4'),
        '2015' = c('H1','H2','H3','H4','H5','H6'),
        '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
        '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
        '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
      ),
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


