# data = oman
# upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay)
# instr_pred_name = 'natural_pred'
# downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay)
# downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02
#
# rain_col_name = 'Rain.Gauge.Measurement'
# gauge_day_type_col_name = 'Gauge.Day.Type'
# upwind_type = 'Upwind'
# downwind_target_type = 'Target'
# downwind_control_type = 'Control'
#
# bootstrap_zero = TRUE

#data is the full data
rain_attr = function(data, upwind_lmm_formula, instr_pred_name, downwind_lmm_formula, downwind_logistic_formula = NULL,
                     rain_col_name,gauge_day_type_col_name, upwind_type, downwind_target_type, downwind_control_type,
                     bootstrap_zero
                     ){
  if(!bootstrap_zero){
    downwind_logistic_formula = NULL
  }

  #Define different subsets of the data
  upwind = (data[,gauge_day_type_col_name] == upwind_type)
  downwind = (data[,gauge_day_type_col_name] %in% c(downwind_target_type, downwind_control_type))
  positive = ( data[,rain_col_name] > 0)

  #Fit Upwind LMM and get the instrumental prediction, then fit Downwind LMM, and Downwind Logistic Model
  fitted_models = fit_upwind_downwind_models(data, upwind_lmm_formula, instr_pred_name, downwind_lmm_formula, downwind_logistic_formula, upwind, downwind, positive)

  #Compute Point Estimates for Attribution - using Ray Winsorize or Proposed Estimates

  #Perform Bootstrap Inference on the attribution estimates (could use parallelization)

  #Perform Permutation Inference on the attribution estimates (could use parallelization)


  return(list(
    fitted_models = list(
      upwind_lmm_fit = upwind_lmm_fit,
      downwind_lmm_fit = downwind_lmm_fit
    )
  ))
}

#data is the full data
fit_upwind_downwind_models = function(data, upwind_lmm_formula, instr_pred_name, downwind_lmm_formula, downwind_logistic_formula = NULL, upwind, downwind, positive){
  upwind_lmm_fit = lme4::lmer(upwind_lmm_formula, data = data[upwind & positive,])
  data = cbind(data, predict(upwind_lmm_fit, data, re.form = NA))
  colnames(data)[ncol(data)] = instr_pred_name
  downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = data[downwind & positive,])
  if(is.null(downwind_logistic_formula)){
    downwind_logistic_fit = NULL
  }else{
    downwind_logistic_fit = glm(downwind_logistic_formula, data = data[downwind,], family = 'binomial')
  }
  return(list(
    upwind_lmm_fit = upwind_lmm_fit,
    downwind_lmm_fit = downwind_lmm_fit,
    downwind_logistic_fit = downwind_logistic_fit
  ))
}
