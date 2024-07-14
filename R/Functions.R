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
# x_downwind_name = c('Gauge.Elevation','natural_pred')
# target_only = TRUE
# attr_type = 'Ray Winsorize'
#
# bootstrap_zero = TRUE

#data is the full data
rain_attr = function(data, upwind_lmm_formula, instr_pred_name, downwind_lmm_formula, downwind_logistic_formula = NULL,
                     rain_col_name,gauge_day_type_col_name, upwind_type, downwind_target_type, downwind_control_type,
                     attr_type, x_downwind_name, target_only,
                     bootstrap_zero
                     ){
  if(!bootstrap_zero){
    downwind_logistic_formula = NULL
  }

  #Define different subsets of the data
  upwind = (data[,gauge_day_type_col_name] == upwind_type)
  downwind = (data[,gauge_day_type_col_name] %in% c(downwind_target_type, downwind_control_type))
  positive = ( data[,rain_col_name] > 0)
  downwind_positive_target = data[downwind & positive, gauge_day_type_col_name] == downwind_target_type
  downwind_positive_control = data[downwind & positive, gauge_day_type_col_name] == downwind_control_type

  #Fit Upwind LMM and get the instrumental prediction, then fit Downwind LMM, and Downwind Logistic Model
  fitted_models = fit_upwind_downwind_models(data, upwind_lmm_formula, instr_pred_name, downwind_lmm_formula, downwind_logistic_formula, upwind, downwind, positive)
  data = fitted_models$data

  #Compute Point Estimates for Attribution - using Ray Winsorize or Proposed Estimates
  hatattr = attr_est(attr_type, data, downwind, positive, rain_col_name, gauge_day_type_col_name, downwind_target_type, downwind_control_type,
                     x_downwind_name, target_only = FALSE, downwind_lmm_fit = fitted_models$downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)


  #Perform Bootstrap Inference on the attribution estimates (could use parallelization) - maybe directly use lme4::lmer() and glm() directly instead of using fit_upwind_downwind_models() in each bootstrap run

  #Perform Permutation Inference on the attribution estimates (could use parallelization) - maybe directly use lme4::lmer() directly instead of using fit_upwind_downwind_models() in each permutation run


  return(list(
    fitted_models = list(
      upwind_lmm_fit = upwind_lmm_fit,
      downwind_lmm_fit = downwind_lmm_fit
    )
  ))
}


#Note that hatu should be for downwind positive (i,t) regardless of target_only
#Optional input arguments: hatalphabeta, hatu
#Note that for attr_type == 'Proposed', the hatSigma_beta matrix is ALWAYS obtained from downwind_lmm_fit
attr_est = function(attr_type, data, downwind, positive, rain_col_name, gauge_day_type_col_name, downwind_target_type, downwind_control_type,
                    x_downwind_name, target_only, downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL){
  if(target_only){
    downwind_positive_useful_row = data[downwind & positive, gauge_day_type_col_name] %in%  c(downwind_target_type)
  }else{
    downwind_positive_useful_row = data[downwind & positive, gauge_day_type_col_name] %in%  c(downwind_target_type,downwind_control_type)
  }

  y_vec = data[downwind & positive,rain_col_name][downwind_positive_useful_row]
  x_z_mat = model.matrix(downwind_lmm_fit, data = data, type = 'fixed')[downwind_positive_useful_row,]

  if(is.null(hatalphabeta)){
    hatalpha_downwind = lme4::fixef(downwind_lmm_fit)[c('(Intercept)',x_downwind_name)]
    hatbeta_downwind = lme4::fixef(downwind_lmm_fit)[setdiff(names(lme4::fixef(downwind_lmm_fit)), c('(Intercept)',x_downwind_name))]
  }else{
    hatalpha_downwind = hatalphabeta[c('(Intercept)',x_downwind_name)]
    hatbeta_downwind = hatalphabeta[setdiff(names(hatalphabeta), c('(Intercept)',x_downwind_name))]
  }



  if(attr_type == 'Ray Winsorize'){
    if(is.null(hatu)){
      hatu = predict(downwind_lmm_fit, newdata = data[downwind & positive,], random.only = TRUE)
    }
    #compute log_hatw = (1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t, for PDR (i,t) or PDR Target (i,t)
    log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row])

    #compute log_hatd = z_it %*% hatbeta, for PDR (i,t) or PDR Target (i,t)
    log_hatd = as.vector(x_z_mat[,setdiff(colnames(x_z_mat),c('(Intercept)',x_downwind_name))] %*% hatbeta_downwind )

    #compute log_haty_naive = (1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t + z_it %*% hatbeta, for PDR (i,t) or PDR Target (i,t)
    log_haty_naive = log_hatw + log_hatd

    #compute mu = (1/n) * sum of  (y_{i't'} / haty_{i't'}^naive ), where the sum is across PDR (i,t) or PDR Target (i,t), and n is the number of PDR (i,t) or PDR Target (i,t)  observation.
    mu = mean( y_vec / exp(log_haty_naive) )

    #compute m = Var{(1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t} / Var{z_it %*% hatbeta}, where the Var is across PDR (i,t) or PDR Target (i,t).
    m = (var(log_hatw))/var(log_hatd)

    #compute lambda based on the formula in Ray's ISR paper
    lambda = 1 + (
      (sqrt( (1 + m)^2 + 4 * (mu - 1) * m  ) - (1 + m)) / (2 * m)
    )

    #compute hatE_it + 1 = min{2, lambda * exp(z_it %*% hatbeta)}, for PDR (i,t) or PDR Target (i,t)
    hatE_plus_one = pmin(2, lambda*exp(log_hatd))

    #hatR_it = y_it / (hatE_it + 1), for PDR (i,t) or PDR Target (i,t)
    hatR =  y_vec / hatE_plus_one

    #compute hatA_it = y_it - hatR_it, for PDR (i,t) or PDR Target (i,t)
    hatA = y_vec - hatR
  }

  if(attr_type == 'Proposed'){
    #compute log_hatd = z_it %*% hatbeta, for PDR (i,t) or PDR Target (i,t)
    log_hatd = as.vector(x_z_mat[,setdiff(colnames(x_z_mat),c('(Intercept)',x_downwind_name))] %*% hatbeta_downwind )

    #Sigma_full is the full covariance matrix of (hatintercept, hatalpha, hatbeta)
    Sigma_full = as.matrix(vcov(downwind_lmm_fit))

    #Sigma_beta = covariance matrix of only hatbeta
    Sigma_beta = Sigma_full[- which(rownames(Sigma_full) %in% c('(Intercept)', x_downwind_name)), -  which(rownames(Sigma_full) %in% c('(Intercept)',x_downwind_name))]

    #Extracting the z_it matrix, for PDR (i,t) or PDR Target (i,t)
    z_mat = x_z_mat[,-which(colnames(x_z_mat) %in%  c('(Intercept)',x_downwind_name) )]

    #Compute the vector of z_it^T %*% Sigma_beta %*% z_it for all PDR (i,t) or PDR Target (i,t)
    zT_Sigma_beta_z = diag(z_mat %*% Sigma_beta %*% t(z_mat))

    lambda = exp(0.5 * zT_Sigma_beta_z)

    #Compute hatR_it = y_it * exp(-z_it %*% beta ) * exp(-(1/2) z_it^T %*% Sigma_beta %*% z_it ), for PDR (i,t) or PDR Target (i,t)
    hatR = y_vec / (exp(log_hatd) * lambda )

    #compute hatA_it = y_it - hatR_it, for PDR (i,t) or PDR Target (i,t)
    hatA = y_vec - hatR
  }

  if(target_only){
    #Ensure that the resulting hatR and hatA has the same length as the number of PDR (i,t), where the values of hatR_it are set to y_it for PDR Control (i,t) and hatA_it are set to 0 for PDR Control (i,t)
    temp = rep(0, sum(downwind & positive))
    temp[downwind_positive_useful_row] = hatA
    hatA = temp
    hatR = data[downwind & positive,rain_col_name]- hatA
  }

  return(list(
    #hatR_it and hatA_it for PDR (i,t), regardless of treatment_only value
    hatR = hatR, hatA = hatA,
    lambda = lambda,
    #attribution quantities, where the sum is either over PDR (i,t) or PDR Target (i,t)
    tor = sum(y_vec),
    tlr = sum(hatR[downwind_positive_useful_row]),
    tac = sum(hatA[downwind_positive_useful_row]),
    apo = sum(hatA[downwind_positive_useful_row])/sum(y_vec),
    apl = sum(hatA[downwind_positive_useful_row])/sum(hatR[downwind_positive_useful_row])
  ))
}

#data is the full data
#optional input arguments: downwind_logistic_formula
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
    downwind_logistic_fit = downwind_logistic_fit,
    data = data
  ))
}
