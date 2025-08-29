# R script to demonstrate functionality of RainAttr package


rm(list=ls())
load('data/oman.rda')
load('data/ionizer_operation.rda')
load('data/gaugeday_downwind.rda')
source('R/Functions.R')


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
  attr_type = 'ThoEtAl',

  #Vector of variable names in downwind_lmm_formula to identify non-ionizer related covariates, whose effects are not included in the calculation of attribution and SATE
  x_downwind_name = c('Year...2013' , 'Year...2014' , 'Year...2016' , 'Year...2017' , 'Year...2018', 'Gauge.Elevation...1km', 'Gauge.Elevation...1km.1', 'natural_pred'),

  #Logical indicating whether the attribution estimates should be based on only 'Target' observations (TRUE), or based on both 'Target' and 'Control' observations (FALSE)
  target_only = FALSE,

  #Indicator for whether to perform bootstrap inference on attribution and SATE
  bootstrap = FALSE,

  #Specification of various option for bootstrap, e.g., bootstrap_option(bootstrap_type = 'PREB1'), currently support bootstrap_type = 'PREB0', 'PREB1', 'PREB2', 'REB0', 'REB1', 'REB2', 'MREB1'
  bootstrap_option = bootstrap_option(),

  #Indicator for whether to perform permutation inference on attribution and SATE
  permutation = FALSE,

  #Specification of various option for permutation,  e.g., whether to permute the operating states between ionizers, between days, between gaugedays
  permutation_option = permutation_option()
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

#### Catching user errors ####

# Inconsistency between instrumental prediction variable name vs. its use in downwind_lmm_formula
# rain_attr(
#   data = oman,
#   upwind_lmm_formula = LogRain ~ Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#   #Variable name for storing instrumental prediction from upwind LMM
#   instr_pred_name = 'natural_pred',
#   instr_pred_type = 'Unconditional',
#   #LMM formula fitted to the downwind observations
#   downwind_lmm_formula = LogRain ~ natural_pred_typo + Year...2013 + Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation...1km:Target.H.01 + Gauge.Elevation...1km:Target.H.02 + + Gauge.Elevation...1km.1:Target.H.01 + Gauge.Elevation...1km.1:Target.H.02 + (1|TrialDay),
#   downwind_logistic_formula = NULL,
#   downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#   rain_col_name = 'Rain.Gauge.Measurement',
#   upwind_subset = Gauge.Day.Type == 'Upwind' & Year!= 2013,
#   downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#   downwind_target_subset = Gauge.Day.Type == 'Target',
#   downwind_control_subset = Gauge.Day.Type == 'Control',
#   positive_subset = Rain.Gauge.Measurement > 0,
#   attr_type = 'ThoEtAl',
#   x_downwind_name = c('Year...2013' , 'Year...2014' , 'Year...2016' , 'Year...2017' , 'Year...2018', 'Gauge.Elevation...1km', 'Gauge.Elevation...1km.1', 'natural_pred'),
#   target_only = FALSE,
#   bootstrap = FALSE,
#   bootstrap_option = bootstrap_option(),
#   permutation = FALSE,
#   permutation_option = permutation_option()
# )

# Invalid attribution type requested by user
# rain_attr(
#   data = oman,
#   upwind_lmm_formula = LogRain ~ Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#   instr_pred_name = 'natural_pred',
#   instr_pred_type = 'Unconditional',
#   downwind_lmm_formula = LogRain ~ Year...2013 + Year...2014 +  Year...2016 + Year...2017 + Year...2018 + Gauge.Elevation...1km + Gauge.Elevation...1km.1 + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation...1km:Target.H.01 + Gauge.Elevation...1km:Target.H.02 + + Gauge.Elevation...1km.1:Target.H.01 + Gauge.Elevation...1km.1:Target.H.02 + (1|TrialDay),
#   downwind_logistic_formula = NULL,
#   downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#   rain_col_name = 'Rain.Gauge.Measurement',
#   upwind_subset = Gauge.Day.Type == 'Upwind' & Year!= 2013,
#   downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#   downwind_target_subset = Gauge.Day.Type == 'Target',
#   downwind_control_subset = Gauge.Day.Type == 'Control',
#   positive_subset = Rain.Gauge.Measurement > 0,
#   #Invalid attr_type requested
#   attr_type = 'Anytype',
#   x_downwind_name = c('Year...2013' , 'Year...2014' , 'Year...2016' , 'Year...2017' , 'Year...2018', 'Gauge.Elevation...1km', 'Gauge.Elevation...1km.1', 'natural_pred'),
#   target_only = FALSE,
#   bootstrap = FALSE,
#   bootstrap_option = bootstrap_option(),
#   permutation = FALSE,
#   permutation_option = permutation_option()
# )

#### Bootstrap Inference ####
set.seed(123)

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
  bootstrap_option = bootstrap_option(
    #Number of bootstrap runs
    B_bootstrap = 3,
    #Type of bootstrap: 'PREB0', 'PREB1', 'PREB2', 'REB0', 'REB1', 'REB2', 'MREB1'
    bootstrap_type = 'PREB1',
    #Whether to bootstrap zeros (i.e., rainfall event indicators)
    bootstrap_zero = F,
    #Probability threshold used for bootstrapping rainfall event indicators
    positive_prob_threshold = NULL,
    #Whether to discretize the bootstrapped rainfall (at raw scale)
    discretize_rain = F,
    #Whether to winsorize individual bootstrapped rainfall to ensure each of them is <= 175mm
    winsorize_individual_rain = F,
    #Whether to winsorize the total bootstrapped rainfall (for all downwind observations) to ensure it is within (6000mm,60000mm)
    winsorize_total_rain = F,
    #Confidence level of bootstrap CI
    CI_level = 0.95
  ),
  permutation = FALSE,
  permutation_option = permutation_option()
)

#Bootstrap CI result for attribution and SATE, as well as the parameters of all fitted models e.g., downwind lmm model, downwind propensity score model
testing_PREB1$bootstrap_CI_result

#Bootstrap p-value result (proportion of bootstrapped estimates that are < 0) for attribution and SATE, as well as the parameters of all fitted models (but less meaningful for these parameters)
testing_PREB1$bootstrap_p_value_result

#Bootstrap distributions for attribution and SATE, as well as the parameters of all fitted models, and the bootstrapped responses (could include offset term) of the downwind_lmm_formula
names(testing_PREB1$bootstrap_result)
testing_PREB1$bootstrap_result[-8]

#Plots of bootstrap distribution for attribution and SATE
ggpubr::ggarrange(plotlist =testing_PREB1$bootstrap_plot_result$hatattr)
ggpubr::ggarrange(plotlist = testing_PREB1$bootstrap_plot_result$hatsate)



#### Permutation Inference ####
set.seed(123)
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
  permutation_option = permutation_option(
    # Number of permutation run
    B_permutation = 5,
    # Whether to permute the elements within each row of ionizer_operation_input
    permute_between_ionizer = T,
    # Whether to permute the rows of ionizer_operation_input
    permute_all_ionizers_between_day = T,
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
    ionizer_operation_day_column_name = 'TrialDay'
  )
)

#Permutation p-value (proportion of permuted estimates that are greater than or equal to original estimate) for attribution and SATE
testing_perm$permutation_p_value_result

#Permutation distribution for attribution and SATE
testing_perm$permutation_result

#Plots of permutation distribution for attribution and SATE
ggpubr::ggarrange(plotlist = testing_perm$permutation_plot_result$hatattr)
ggpubr::ggarrange(plotlist = testing_perm$permutation_plot_result$hatsate)
