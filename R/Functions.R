load("data/oman.rda")
data = oman
upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay)
instr_pred_name = 'natural_pred'
downwind_lmm_formula = LogRain - natural_pred ~ Gauge.Elevation  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay)
downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02

rain_col_name = 'Rain.Gauge.Measurement'

#Old way
# gauge_day_type_col_name = 'Gauge.Day.Type'
# upwind_type = 'Upwind'
# downwind_target_type = 'Target'
# downwind_control_type = 'Control'


##New way
upwind_subset = Gauge.Day.Type == 'Upwind'
downwind_subset = Gauge.Day.Type  %in% c('Target','Control')
downwind_target_subset = Gauge.Day.Type == 'Target'
downwind_control_subset = Gauge.Day.Type == 'Control'

x_downwind_name = c('Gauge.Elevation','natural_pred')
target_only = FALSE
attr_type = 'Ray Winsorize'

bootstrap_zero = TRUE



#data is the full data
#To use instr_pred or natural_pred as an offset, just specify the response to be LogRain - natural_pred in downwind_lmm_formula
#To fit first-stage model to Control^#, jut specify upwind_subset = Gauge.Day.Type == 'Control'

rain_attr = function(data, upwind_lmm_formula, instr_pred_name, instr_pred_type,
                     downwind_lmm_formula, downwind_logistic_formula = NULL, downwind_propensity_formula,
                     rain_col_name,
                     upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset,
                     attr_type, x_downwind_name, target_only,
                     bootstrap_zero
                     ){


  if(!bootstrap_zero){
    downwind_logistic_formula = NULL
  }

  if(!instr_pred_name %in% all.vars(downwind_lmm_formula)){
    stop("instr_pred_name cannot be found in downwind_lmm_formula")
  }

  if(formula.tools::lhs(as.formula(gsub("[()]", "", downwind_propensity_formula))) != substitute(downwind_target_subset)){
    stop("The definition of Target provided in  downwind_target_subset is not consistent with the LHS of downwind_propensity_formula")
  }

  if(mean(x_downwind_name %in% all.vars(formula.tools::rhs(downwind_lmm_formula)) )!=1 ){
    stop("At least one variable in x_downwind_name cannot be found on RHS of downwind_lmm_formula")
  }


  #Define different subsets of the data
  #Binary indicator of length N, indicating whether or not each observation is an upwind observation
  upwind = eval(substitute(upwind_subset), data, parent.frame())

  #Binary indicator of length N, indicating whether or not each observation is a downwind observation
  downwind = eval(substitute(downwind_subset), data, parent.frame())

  #Binary indicator of length N, indicating whether or not each observation has positive rainfall
  #Could consider to replace this by positive = (!is.na(data[,all.vars(downwind_lmm_formula)[1]])), which allow us to drop rain_col_name, but we still need rain_col_name to compute the attribution estimate anyway
  positive = ( data[,rain_col_name] > 0)

  #Binary indicator of length N_downwind_positive, indicating whether or not each downwind positive observation is target observation
  downwind_positive_target = eval(substitute(downwind_target_subset), data[downwind & positive,], parent.frame())

  #Binary indicator of length N_downwind_positive, indicating whether or not each downwind positive observation is control observation
  downwind_positive_control = eval(substitute(downwind_control_subset), data[downwind & positive,], parent.frame())

  #Fit Upwind LMM and get the instrumental prediction, then fit Downwind LMM, and Downwind Logistic Model
  fitted_models = fit_upwind_downwind_models(data, upwind_lmm_formula, instr_pred_name, instr_pred_type, downwind_lmm_formula, downwind_logistic_formula, upwind, downwind, positive)
  data = fitted_models$data

  #Compute Point Estimates for Attribution - using Ray Winsorize or Proposed Estimates
  hatattr = attr_est(attr_type, data, downwind, positive, rain_col_name, downwind_positive_target, downwind_positive_control,
                     x_downwind_name, target_only = FALSE, downwind_lmm_fit = fitted_models$downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)




  #Compute Points Estimates for SATE - using different types of SATE estimates
  z_downwind_name = setdiff(names(lme4::fixef(fitted_models$downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
  downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)
  hatsate = sate_est(data, downwind, positive, rain_col_name, downwind_positive_target, downwind_positive_control, downwind_propensity_formula, downwind_separate_formula,
                     x_downwind_name, downwind_lmm_fit = fitted_models$downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

  #Perform Bootstrap Inference on the attribution estimates (could use parallelization) - maybe directly use lme4::lmer() and glm() directly instead of using fit_upwind_downwind_models() in each bootstrap run
  #Bootstrap function can follow similar attribute_bootstrap() in D:\Postdoc\Simulation\Replicate ISR Results\Bootstrap Analysis with generate_zero_T and scaled_h_sampling and Correct Scaling REB1 using Oman Data.R
  #as well as D:\Postdoc\Bootstrap Paper\R Codes\Functions_realdata.R
  #Also look at Overleaf/Rainfall Enhancement/Ray's implementation.tex

  #Bootstrap function should allow for the choice of bootstrap_type (REB0/1/2, PREB0/1/2, MREB-1), as well as whether or not to bootstrap zero.
  #When bootstrap_zero = T, need to check Ray's original code and my implementation of MQ bootstrap to see how we get bootstrap distribution of SATE - more specifically, did we refit the propensity logistic model for each bootstrap dataset?

  #For bootstrapping, need to be careful when we are using instr_pred as offset - the generation of y_b data should be different in this case, if not, maybe need to keep in mind we are modelling LogRain - natural_pred


  #Perform Permutation Inference on the attribution estimates (could use parallelization) - maybe directly use lme4::lmer() directly instead of using fit_upwind_downwind_models() in each permutation run


  return(list(
    fitted_models = list(
      upwind_lmm_fit = fitted_models$upwind_lmm_fit,
      downwind_lmm_fit = fitted_models$downwind_lmm_fit,
      downwind_logistic_fit = fitted_models$downwind_logistic_fit
    ),
    hatattr = hatattr,
    hatsate = hatsate
  ))
}



#Note that hatu should be for downwind positive (i,t) regardless of target_only
#Optional input arguments: hatalphabeta, hatu
#Note that for attr_type == 'Proposed', the hatSigma_beta matrix is ALWAYS obtained from downwind_lmm_fit
#When we consider hatbeta from MQ, can add another option to use either hatbeta_{0.5} or hatbeta_{conditional}
attr_est = function(attr_type, data, downwind, positive, rain_col_name, downwind_positive_target, downwind_positive_control,
                    x_downwind_name, target_only, downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL){
  if(target_only){
    downwind_positive_useful_row = downwind_positive_target
  }else{
    downwind_positive_useful_row = (downwind_positive_target | downwind_positive_control)
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

#CONTINUE FROM HERE: check the expressions of different sate estimates below, and cross-check their estimated values with my previous results (either bootstrap paper or the MQ linear regression) or Ray's previous results
#It seems like our sate estimates here are different from those using MQ linear regression, this is because those settings in MQ Linear regression only considers LMM in the first stage but the second stage is ALWAYS MQ!

#Compute different types of SATE estimate in Chambers et al. (2022)
#When we consider hatbeta from MQ, can add another option to use either hatbeta_{0.5} or hatbeta_{conditional}
#Note that when the instr_pred or natural_pred is used as an offset term, the sate.ipw is computed using LogRain - natural_pred, instead of LogRain only
sate_est = function(data, downwind, positive, rain_col_name, downwind_positive_target, downwind_positive_control, downwind_propensity_formula, downwind_separate_formula,
                    x_downwind_name, downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL){
  downwind_propensity_fit = glm(downwind_propensity_formula, family = binomial, data = data[downwind & positive,])
  hatpi = predict(downwind_propensity_fit, type = "response")
  hatw_1 = (1/hatpi)/( sum( (1/hatpi) * as.numeric(downwind_propensity_fit$y) )   )
  hatw_0 = (1/(1-hatpi))/( sum( (1/ (1-hatpi)  ) * (1-as.numeric(downwind_propensity_fit$y))  )   )


  x_z_mat = model.matrix(downwind_lmm_fit, data = data, type = 'fixed')
  z_mat = x_z_mat[,-which(colnames(x_z_mat) %in%  c('(Intercept)',x_downwind_name) )]

  if(is.null(hatalphabeta)){
    hatalpha_downwind = lme4::fixef(downwind_lmm_fit)[c('(Intercept)',x_downwind_name)]
    hatbeta_downwind = lme4::fixef(downwind_lmm_fit)[setdiff(names(lme4::fixef(downwind_lmm_fit)), c('(Intercept)',x_downwind_name))]
  }else{
    hatalpha_downwind = hatalphabeta[c('(Intercept)',x_downwind_name)]
    hatbeta_downwind = hatalphabeta[setdiff(names(hatalphabeta), c('(Intercept)',x_downwind_name))]
  }


  sate.mb = sum( as.numeric(downwind_propensity_fit$y) * as.vector(z_mat %*% hatbeta_downwind))/sum( as.numeric(downwind_propensity_fit$y)  )
  sate.ipw = sum( hatw_1 * as.numeric(downwind_propensity_fit$y) * lme4::getME(downwind_lmm_fit, 'y')  ) - sum( hatw_0 * ( 1- as.numeric(downwind_propensity_fit$y)) * lme4::getME(downwind_lmm_fit, 'y'  ) )
  sate.ipw.l = sum( hatw_1 * as.numeric(downwind_propensity_fit$y) * as.vector(z_mat %*% hatbeta_downwind)  )


  #Fit separate LMM to downwind_positive_target and downwind_positive_control
  downwind_positive_target_lmm_fit = lme4::lmer(downwind_separate_formula, data = data[downwind & positive,][downwind_positive_target,])
  downwind_positive_control_lmm_fit = lme4::lmer(downwind_separate_formula, data = data[downwind & positive,][downwind_positive_control,])
  hatm_1 = predict(downwind_positive_target_lmm_fit, newdata = data[downwind & positive,], re.form = NA)
  hatm_0 = predict(downwind_positive_control_lmm_fit, newdata = data[downwind & positive,], re.form = NA)
  sate.ipw.ma = mean(hatm_1) - mean(hatm_0) + sum(hatw_1 * as.numeric(downwind_propensity_fit$y) * (lme4::getME(downwind_lmm_fit, 'y') - hatm_1)  ) - sum( hatw_0 * (1 - as.numeric(downwind_propensity_fit$y) ) * (lme4::getME(downwind_lmm_fit, 'y'  ) - hatm_0)  )
  sate.aipw = sum( hatw_1 * ( as.numeric(downwind_propensity_fit$y) * lme4::getME(downwind_lmm_fit, 'y') - ( (as.numeric(downwind_propensity_fit$y) - hatpi ) * hatm_1   )  )  ) - sum( hatw_0 * ( (1 - as.numeric(downwind_propensity_fit$y)) *  lme4::getME(downwind_lmm_fit, 'y'  ) -  ( (as.numeric(downwind_propensity_fit$y) - hatpi ) * hatm_0   )   )   )

  #CONTINUE FROM HERE: check why we can't recover the relationship between AIPW and IPW.MA as claimed below equation (7) of Ray causal paper
  check = sate.aipw - sate.ipw.ma
  check2 = sum( hatm_1 * ( 1/(sum( (1/hatpi) *as.numeric(downwind_propensity_fit$y)  ) ) - 1/sum(downwind & positive)   ) ) -  sum(hatm_0 * ( 1/( sum( (1/(1-hatpi)) * (1 - as.numeric(downwind_propensity_fit$y))  ) ) - 1/sum(downwind & positive)  ) )
  browser()

  return(list(
    sate.mb = sate.mb,
    sate.ipw = sate.ipw,
    sate.ipw.l = sate.ipw.l,
    sate.ipw.ma = sate.ipw.ma,
    sate.aipw = sate.aipw
  ))
}

#data is the full data
#optional input arguments: downwind_logistic_formula
fit_upwind_downwind_models = function(data, upwind_lmm_formula, instr_pred_name, instr_pred_type, downwind_lmm_formula, downwind_logistic_formula = NULL, upwind, downwind, positive){
  upwind_lmm_fit = lme4::lmer(upwind_lmm_formula, data = data[upwind & positive,])

  if(!instr_pred_type %in% c('Unconditional','Conditional')){
    stop("instr_pred_type should be either 'Unconditional' or 'Conditional'")
  }

  if(instr_pred_type == 'Unconditional'){
    data = cbind(data, predict(upwind_lmm_fit, data, re.form = NA))
  }
  if(instr_pred_type == 'Conditional'){
    data = cbind(data, predict(upwind_lmm_fit, data, re.form = NULL, allow.new.levels = T))
  }
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

# Note that this function requires downwind_lmm_formula to not use the interaction syntax such as x1*x2, but instead it should always use x1 + x2 + x1:x2
#Function to remove some variables from a given formula - this is used to obtain the formula for fitting separate LMM to downwind_positive_target and downwind_positive_control
remove_fixed_terms <- function(input_formula, vars_to_remove){
  # Flatten input_formula to a single string
  rhs_str = paste(deparse(formula.tools::rhs(input_formula)), collapse = "")
  rhs_terms = trimws(unlist(strsplit(rhs_str,'\\+')))
  rhs_terms = rhs_terms[!rhs_terms %in% vars_to_remove]

  new_rhs_str = paste(rhs_terms, collapse = ' + ')
  lhs_str = paste(deparse(formula.tools::lhs(input_formula)), collapse = "")

  # Build new input_formula
  return(as.formula(paste(lhs_str, "~", new_rhs_str)))
}

#For checking: apo should be 0.111206, apl should be 0.1251201 for attr_type = 'Ray Winsorize'
# apo  = 0.06265952, apl =  0.0668482 for attr_type = 'Proposed'
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


data = oman
upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay)
instr_pred_name = 'natural_pred'
downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay)
downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02
downwind_propensity_formula = Gauge.Day.Type == 'Target' ~ Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure
rain_col_name = 'Rain.Gauge.Measurement'
#

asd  = rain_attr(data = oman,
                 upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
                 instr_pred_name = 'natural_pred',
                 instr_pred_type = 'Unconditional',
                 downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
                 downwind_logistic_formula = NULL,
                 downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
                 rain_col_name = 'Rain.Gauge.Measurement',
                 upwind_subset = Gauge.Day.Type == 'Upwind',
                 downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
                 downwind_target_subset = Gauge.Day.Type == 'Target',
                 downwind_control_subset = Gauge.Day.Type == 'Control',
                 attr_type = 'Ray Winsorize',
                 x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
                 target_only = FALSE,
                 bootstrap_zero = FALSE)

asd$hatsate$sate.ipw
