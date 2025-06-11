# load("data/oman.rda")
# data = oman
# upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay)
# instr_pred_name = 'natural_pred'
# downwind_lmm_formula = LogRain - natural_pred ~ Gauge.Elevation  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay)
# downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02
#
# rain_col_name = 'Rain.Gauge.Measurement'
#
# #Old way
# # gauge_day_type_col_name = 'Gauge.Day.Type'
# # upwind_type = 'Upwind'
# # downwind_target_type = 'Target'
# # downwind_control_type = 'Control'
#
#
# ##New way
# upwind_subset = Gauge.Day.Type == 'Upwind'
# downwind_subset = Gauge.Day.Type  %in% c('Target','Control')
# downwind_target_subset = Gauge.Day.Type == 'Target'
# downwind_control_subset = Gauge.Day.Type == 'Control'
#
# x_downwind_name = c('Gauge.Elevation','natural_pred')
# target_only = FALSE
# attr_type = 'Ray Winsorize'
#
# bootstrap_zero = TRUE



#data is the full data
#To use instr_pred or natural_pred as an offset, just specify the response to be LogRain - natural_pred in downwind_lmm_formula
#To fit first-stage model to Control^#, jut specify upwind_subset = Gauge.Day.Type == 'Control'
rm(list=ls())
load('data/oman.rda')
load('data/ionizer_operation.rda')
load('data/gaugeday_downwind.rda')
rain_attr = function(data, upwind_lmm_formula, instr_pred_name, instr_pred_type,
                     downwind_lmm_formula, downwind_logistic_formula = NULL, downwind_propensity_formula,
                     rain_col_name,
                     upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset,
                     attr_type, x_downwind_name, target_only,
                     bootstrap = FALSE, bootstrap_option = NULL,
                     permutation = FALSE, permutation_option = NULL
                     ){
  #TODO: Add checks for permutation_option such as checking if permutation_option$ionizer_operation_year_column_name and permutation_option$ionizer_operation_day_column_name can be found in the colnames(ionizer_operation)

  if(!instr_pred_name %in% all.vars(downwind_lmm_formula)){
    stop("instr_pred_name cannot be found in downwind_lmm_formula")
  }

  if(formula.tools::lhs(as.formula(gsub("[()]", "", downwind_propensity_formula))) != substitute(downwind_target_subset)){
    stop("The definition of Target provided in  downwind_target_subset is not consistent with the LHS of downwind_propensity_formula")
  }

  if(mean(x_downwind_name %in% all.vars(formula.tools::rhs(downwind_lmm_formula)) )!=1 ){
    stop("At least one variable in x_downwind_name cannot be found on RHS of downwind_lmm_formula")
  }

  if(!attr_type %in% c('Ray Winsorize','Proposed', 'No')){
    stop("attr_type should be one of 'Ray Winsorize','Proposed', or 'No''")
  }

  if(permutation){
    if(!permutation_option$ionizer_operation_year_column_name %in% colnames(permutation_option$ionizer_operation) ){
      stop("The columns of permutation_option$ionizer_operation do not contain permutation_option$ionizer_operation_year_column_name")
    }

    if(!permutation_option$ionizer_operation_day_column_name %in% colnames(permutation_option$ionizer_operation) ){
      stop("The columns of permutation_option$ionizer_operation do not contain permutation_option$ionizer_operation_day_column_name")
    }

    if(!permutation_option$ionizer_operation_day_column_name %in% colnames(data) ){
      stop("The columns of data do not contain permutation_option$ionizer_operation_day_column_name")
    }

    if(mean(permutation_option$data_target_column_names %in% colnames(data)) != 1){
      stop("At least one element of permutation_option$data_target_column_names is cannot be found in the column names of data")
    }

    if(length(unique(permutation_option$ionizer_operation[,permutation_option$ionizer_operation_day_column_name])) != length(permutation_option$ionizer_operation[,permutation_option$ionizer_operation_day_column_name]) ){
      stop('permutation_option$ionizer_operation[,permutation_option$ionizer_operation_day_column_name] contains duplicated values')
    }

    if(mean(data[,permutation_option$ionizer_operation_day_column_name] %in% permutation_option$ionizer_operation[,permutation_option$ionizer_operation_day_column_name] ) !=1 ){
      stop('At least one day in data cannot be found in permutation_option$ionizer_operation')
    }

    if(nrow(data)!= nrow(permutation_option$gaugeday_downwind)){
      stop('Number of rows in data does not equal number of rows in permutation_option$gaugeday_downwind')
    }

    if( length(permutation_option$year_ionizer_list) != length(unique(names(permutation_option$year_ionizer_list)))  ){
      stop('permutation_option$year_ionizer_list has duplicated names')
    }

    if( mean( permutation_option$ionizer_operation[,permutation_option$ionizer_operation_year_column_name] %in% names(permutation_option$year_ionizer_list) )!=1 ){
      stop('At least one element of permutation_option$ionizer_operation[,permutation_option$ionizer_operation_year_column_name] cannot be found in the element names of permutation_option$year_ionizer_list')
    }
  }

  #Define different subsets of the data
  #Binary indicator of length N, indicating whether or not each observation is an upwind observation
  upwind_expr = rlang::enquo(upwind_subset)
  upwind = rlang::eval_tidy(upwind_expr, data = data)

  #Binary indicator of length N, indicating whether or not each observation is a downwind observation
  downwind_expr = rlang::enquo(downwind_subset)
  downwind = rlang::eval_tidy(downwind_expr, data = data)

  #Binary indicator of length N, indicating whether or not each observation has positive rainfall
  #Could consider to replace this by positive = (!is.na(data[,all.vars(downwind_lmm_formula)[1]])), which allow us to drop rain_col_name, but we still need rain_col_name to compute the attribution estimate anyway
  positive = ( data[,rain_col_name] > 0)

  #Binary indicator of length N_downwind_positive, indicating whether or not each downwind positive observation is target observation
  downwind_positive_target_expr = rlang::enquo(downwind_target_subset)
  downwind_positive_target = rlang::eval_tidy(downwind_positive_target_expr, data = data[downwind & positive,])

  #Binary indicator of length N_downwind_positive, indicating whether or not each downwind positive observation is control observation
  downwind_positive_control_expr = rlang::enquo(downwind_control_subset)
  downwind_positive_control = rlang::eval_tidy(downwind_positive_control_expr, data = data[downwind & positive,])

  #Fit Upwind LMM and get the instrumental prediction, then fit Downwind LMM, and Downwind Logistic Model
  fitted_models = fit_upwind_downwind_models(data, upwind_lmm_formula, instr_pred_name, instr_pred_type, downwind_lmm_formula, downwind_logistic_formula, upwind, downwind, positive)

  downwind_positive_data = fitted_models$data[downwind & positive, ]
  #Compute Point Estimates for Attribution - using Ray Winsorize or Proposed Estimates
  hatattr = attr_est(attr_type, downwind_positive_data, rain_col_name, downwind_positive_target, downwind_positive_control,
                     x_downwind_name, target_only = target_only, downwind_lmm_fit = fitted_models$downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)




  #Compute Points Estimates for SATE - using different types of SATE estimates
  z_downwind_name = setdiff(names(lme4::fixef(fitted_models$downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
  downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)
  hatsate = sate_est(downwind_positive_data, downwind_positive_target, downwind_positive_control, downwind_propensity_formula, downwind_separate_formula,
                     x_downwind_name, downwind_lmm_fit = fitted_models$downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

  #Perform Bootstrap Inference on the attribution estimates (could use parallelization) - maybe directly use lme4::lmer() and glm() directly instead of using fit_upwind_downwind_models() in each bootstrap run
  #Bootstrap function can follow similar attribute_bootstrap() in D:\Postdoc\Simulation\Replicate ISR Results\Bootstrap Analysis with generate_zero_T and scaled_h_sampling and Correct Scaling REB1 using Oman Data.R
  #as well as D:\Postdoc\Bootstrap Paper\R Codes\Functions_realdata.R
  #Also look at Overleaf/Rainfall Enhancement/Ray's implementation.tex

  #Bootstrap function should allow for the choice of bootstrap_type (REB0/1/2, PREB0/1/2, MREB-1), as well as whether or not to bootstrap zero.
  #When bootstrap_zero = T, need to check Ray's original code and my implementation of MQ bootstrap to see how we get bootstrap distribution of SATE - more specifically, did we refit the propensity logistic model for each bootstrap dataset?

  all_fitted_models = list(
    upwind_lmm_fit = fitted_models$upwind_lmm_fit,
    downwind_lmm_fit = fitted_models$downwind_lmm_fit,
    downwind_logistic_fit = fitted_models$downwind_logistic_fit,
    downwind_propensity_fit = hatsate$fitted_models$downwind_propensity_fit,
    downwind_positive_target_lmm_fit = hatsate$fitted_models$downwind_positive_target_lmm_fit,
    downwind_positive_control_lmm_fit = hatsate$fitted_models$downwind_positive_control_lmm_fit
  )

  #For bootstrapping, need to be careful when we are using instr_pred as offset - the generation of y_b data should be different in this case, if not, maybe need to keep in mind we are modelling LogRain - natural_pred
  if(bootstrap){
    bootstrap_result = bootstrap_downwind(B_bootstrap = bootstrap_option$B_bootstrap,
                                          bootstrap_type = bootstrap_option$bootstrap_type,
                                          bootstrap_zero = bootstrap_option$bootstrap_zero,
                                          positive_prob_threshold = bootstrap_option$positive_prob_threshold,
                                          discretize_rain = bootstrap_option$discretize_rain,
                                          winsorize_individual_rain = bootstrap_option$winsorize_individual_rain,
                                          winsorize_total_rain = bootstrap_option$winsorize_total_rain,
                                          ori_data = fitted_models$data,
                                          downwind = downwind,
                                          ori_positive = positive,
                                          rain_col_name = rain_col_name,
                                          downwind_target_subset = !!downwind_positive_target_expr,
                                          downwind_control_subset = !!downwind_positive_control_expr,
                                          ori_fitted_models = all_fitted_models,
                                          downwind_lmm_formula = downwind_lmm_formula, attr_type = attr_type, x_downwind_name = x_downwind_name, target_only = target_only,
                                          downwind_propensity_formula = downwind_propensity_formula,
                                          ori_attr_est = c(hatattr$apo,hatattr$apl),
                                          ori_sate_est = c(hatsate$estimates$sate.mb, hatsate$estimates$sate.ipw, hatsate$estimates$sate.ipw.l, hatsate$estimates$sate.ipw.ma, hatsate$estimates$sate.aipw)
                                          )

    bootstrap_CI_result = lapply(bootstrap_result[-length(bootstrap_result)], function(x){bootstrap_CI(x,level = bootstrap_option$CI_level)})
    bootstrap_p_value_result = lapply(bootstrap_result[-length(bootstrap_result)], function(x){bootstrap_p_value(x)})
    bootstrap_plot_result = list(
      hatattr = bootstrap_plot(bootstrap_result$hatattr, ori_est = c(hatattr$apo,hatattr$apl)),
      hatsate = bootstrap_plot(bootstrap_result$hatsate, ori_est = c(hatsate$estimates$sate.mb, hatsate$estimates$sate.ipw, hatsate$estimates$sate.ipw.l, hatsate$estimates$sate.ipw.ma, hatsate$estimates$sate.aipw))
    )


  }else{
    bootstrap_result = NULL
    bootstrap_CI_result = NULL
    bootstrap_p_value_result = NULL
    bootstrap_plot_result = NULL
  }




  #Perform Permutation Inference on the attribution estimates (could use parallelization) - maybe directly use lme4::lmer() directly instead of using fit_upwind_downwind_models() in each permutation run

  if(permutation){
    permutation_result = permutation_ionizer(B_permutation = permutation_option$B_permutation,
                                             permute_between_ionizer = permutation_option$permute_between_ionizer,
                                             permute_all_ionizers_between_day = permutation_option$permute_all_ionizers_between_day,
                                             permute_between_gaugeday = permutation_option$permute_between_gaugeday,
                                             ionizer_operation = permutation_option$ionizer_operation,
                                             gaugeday_downwind = permutation_option$gaugeday_downwind,
                                             year_ionizer_list = permutation_option$year_ionizer_list,
                                             data_target_column_names = permutation_option$data_target_column_names,
                                             ionizer_operation_year_column_name = permutation_option$ionizer_operation_year_column_name,
                                             ionizer_operation_day_column_name = permutation_option$ionizer_operation_day_column_name,
                                             data = fitted_models$data,
                                             downwind_lmm_formula = downwind_lmm_formula,
                                             downwind_propensity_formula = downwind_propensity_formula,
                                             attr_type = attr_type,
                                             x_downwind_name = x_downwind_name,
                                             target_only = target_only,
                                             rain_col_name = rain_col_name)
    permutation_p_value_result = list(
      hatattr = permutation_p_value(permutation_result$hatattr, ori_est = c(hatattr$apo,hatattr$apl)),
      hatsate = permutation_p_value(permutation_result$hatsate, ori_est = c(hatsate$estimates$sate.mb, hatsate$estimates$sate.ipw, hatsate$estimates$sate.ipw.l, hatsate$estimates$sate.ipw.ma, hatsate$estimates$sate.aipw))
    )

    permutation_plot_result = list(
      hatattr = permutation_plot(permutation_result$hatattr, ori_est = c(hatattr$apo,hatattr$apl)),
      hatsate = permutation_plot(permutation_result$hatsate, ori_est = c(hatsate$estimates$sate.mb, hatsate$estimates$sate.ipw, hatsate$estimates$sate.ipw.l, hatsate$estimates$sate.ipw.ma, hatsate$estimates$sate.aipw))
    )
  }else{
    permutation_result = NULL
    permutation_p_value_result = NULL
    permutation_plot_result = NULL
  }

  return(list(
    all_fitted_models = all_fitted_models,
    hatattr = hatattr,
    hatsate = hatsate$estimates,
    bootstrap_result = bootstrap_result,
    bootstrap_CI_result = bootstrap_CI_result,
    bootstrap_p_value_result = bootstrap_p_value_result,
    bootstrap_plot_result = bootstrap_plot_result,
    permutation_result = permutation_result,
    permutation_p_value_result = permutation_p_value_result,
    permutation_plot_result = permutation_plot_result,
    #temporary - will be removed later sicne we dont need to return data. Currently included for debugging purposes
    data = fitted_models$data
  ))
}



#Note that hatu should be for downwind positive (i,t) regardless of target_only
#Optional input arguments: hatalphabeta, hatu
#Note that for attr_type == 'Proposed', the hatSigma_beta matrix is ALWAYS obtained from downwind_lmm_fit
#When we consider hatbeta from MQ, can add another option to use either hatbeta_{0.5} or hatbeta_{conditional}
attr_est = function(attr_type, downwind_positive_data, rain_col_name, downwind_positive_target, downwind_positive_control,
                    x_downwind_name, target_only, downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL){
  if(target_only){
    downwind_positive_useful_row = downwind_positive_target
  }else{
    downwind_positive_useful_row = (downwind_positive_target | downwind_positive_control)
  }

  y_vec = downwind_positive_data[downwind_positive_useful_row,rain_col_name]
  x_z_mat = model.matrix(downwind_lmm_fit, data = downwind_positive_data, type = 'fixed')[downwind_positive_useful_row,]

  if(is.null(hatalphabeta)){
    hatalpha_downwind = lme4::fixef(downwind_lmm_fit)[c('(Intercept)',x_downwind_name)]
    hatbeta_downwind = lme4::fixef(downwind_lmm_fit)[setdiff(names(lme4::fixef(downwind_lmm_fit)), c('(Intercept)',x_downwind_name))]
  }else{
    hatalpha_downwind = hatalphabeta[c('(Intercept)',x_downwind_name)]
    hatbeta_downwind = hatalphabeta[setdiff(names(hatalphabeta), c('(Intercept)',x_downwind_name))]
  }



  if(attr_type == 'Ray Winsorize'){
    if(is.null(hatu)){
      hatu = predict(downwind_lmm_fit, newdata = downwind_positive_data, random.only = TRUE)
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

  if(attr_type == 'No'){
    #compute log_hatd = z_it %*% hatbeta, for PDR (i,t) or PDR Target (i,t)
    log_hatd = as.vector(x_z_mat[,setdiff(colnames(x_z_mat),c('(Intercept)',x_downwind_name))] %*% hatbeta_downwind )

    lambda = 1

    #Compute hatR_it = y_it * exp(-z_it %*% beta ) * exp(-(1/2) z_it^T %*% Sigma_beta %*% z_it ), for PDR (i,t) or PDR Target (i,t)
    hatR = y_vec / (exp(log_hatd) * lambda )

    #compute hatA_it = y_it - hatR_it, for PDR (i,t) or PDR Target (i,t)
    hatA = y_vec - hatR
  }

  if(target_only){
    #Ensure that the resulting hatR and hatA has the same length as the number of PDR (i,t), where the values of hatR_it are set to y_it for PDR Control (i,t) and hatA_it are set to 0 for PDR Control (i,t)
    temp = rep(0, nrow(downwind_positive_data))
    temp[downwind_positive_useful_row] = hatA
    hatA = temp
    hatR = downwind_positive_data[,rain_col_name]- hatA
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


#DONE: verified the sate_est() below, by having the same exact SATE estimates as 'Different Variations of 2SLMM on Log Rainfall.R', as well as same reported SATE estimates in Table 6 of JRSSA paper corresponding to LogRain as response.

#Compute different types of SATE estimate in Chambers et al. (2022)
#When we consider hatbeta from MQ, can add another option to use either hatbeta_{0.5} or hatbeta_{conditional}
#Note that when the instr_pred or natural_pred is used as an offset term, the sate.ipw is computed using LogRain - natural_pred, instead of LogRain only
sate_est = function(downwind_positive_data, downwind_positive_target, downwind_positive_control, downwind_propensity_formula, downwind_separate_formula,
                    x_downwind_name, downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL){

  downwind_propensity_fit = glm(downwind_propensity_formula, family = binomial, data = downwind_positive_data)
  hatpi = predict(downwind_propensity_fit, type = "response")

  hatw_1 = (1/hatpi)/( sum( (1/hatpi) * as.numeric(downwind_propensity_fit$y) )   )
  hatw_0 = (1/(1-hatpi))/( sum( (1/ (1-hatpi)  ) * (1-as.numeric(downwind_propensity_fit$y))  )   )


  x_z_mat = model.matrix(downwind_lmm_fit, data = downwind_positive_data, type = 'fixed')
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

  # #Checking relationship between IPW-L and IPW in eq(4) of JRSSA - the expression is correct
  # r.vec = lme4::getME(downwind_lmm_fit, 'y')  - as.vector(z_mat %*% hatbeta_downwind)
  # sate.ipw.l.check1 = sate.ipw - (  sum(hatw_1 *as.numeric(downwind_propensity_fit$y) * r.vec ) - sum(hatw_0 * (1 - as.numeric(downwind_propensity_fit$y)) * r.vec ) )
  # sate.ipw.l - sate.ipw.l.check1
  #
  # r.vec2 = lme4::getME(downwind_lmm_fit, 'y')  - as.vector(z_mat %*% hatbeta_downwind) *as.numeric(downwind_propensity_fit$y)
  # sate.ipw.l.check2 = sate.ipw - (  sum(hatw_1 *as.numeric(downwind_propensity_fit$y) * r.vec2 ) - sum(hatw_0 * (1 - as.numeric(downwind_propensity_fit$y)) * r.vec2 ) )
  # sate.ipw.l - sate.ipw.l.check2
  # browser()

  #Fit separate LMM to downwind_positive_target and downwind_positive_control
  downwind_positive_target_lmm_fit = lme4::lmer(downwind_separate_formula, data = downwind_positive_data[downwind_positive_target,])
  downwind_positive_control_lmm_fit = lme4::lmer(downwind_separate_formula, data = downwind_positive_data[downwind_positive_control,])
  hatm_1 = predict(downwind_positive_target_lmm_fit, newdata = downwind_positive_data, re.form = NA)
  hatm_0 = predict(downwind_positive_control_lmm_fit, newdata = downwind_positive_data, re.form = NA)
  sate.ipw.ma = mean(hatm_1) - mean(hatm_0) + sum(hatw_1 * as.numeric(downwind_propensity_fit$y) * (lme4::getME(downwind_lmm_fit, 'y') - hatm_1)  ) - sum( hatw_0 * (1 - as.numeric(downwind_propensity_fit$y) ) * (lme4::getME(downwind_lmm_fit, 'y'  ) - hatm_0)  )

  #TODO: Need to check when this function is used to compute estimated sate.aipw within each bootstrap run and we are using natural_pred as offset term, should we still follow the equation (7) in JRSSA paper,
  #where we replace y_i with (LogRain_i - natural_pred_i) which is captured by lme4::getME(downwind_lmm_fit, 'y') ) below that returns the response of b_downwind_lmm_fit i.e., bootstrapped_y
  sate.aipw = sum( hatw_1 * ( ( as.numeric(downwind_propensity_fit$y) * lme4::getME(downwind_lmm_fit, 'y') ) - ( (as.numeric(downwind_propensity_fit$y) - hatpi ) * hatm_1   )  )  ) - sum( hatw_0 * ( ( (1 - as.numeric(downwind_propensity_fit$y)) *  lme4::getME(downwind_lmm_fit, 'y'  )  ) -  ( (as.numeric(downwind_propensity_fit$y) - hatpi ) * hatm_0   )   )   )

  # #Checking the equation below eq(7) of JRSSA
  # check = sate.aipw - sate.ipw.ma
  # #This expression provided in the equation below eq(7) of JRSSA paper is probably incorrect
  check2 = sum( hatm_1 * ( 1/(sum( (1/hatpi) *as.numeric(downwind_propensity_fit$y)  ) ) - 1/ nrow(downwind_positive_data)   ) ) -  sum(hatm_0 * ( 1/( sum( (1/(1-hatpi)) * (1 - as.numeric(downwind_propensity_fit$y) )  ) ) - 1/ nrow(downwind_positive_data)  ) )
  # #This expression is the correct expression - see its derivation in iPad's 'Relationship between AIPW and IPW-MA'
  # check3 = sum( hatm_1 * ( 1/(sum( (1/hatpi) *as.numeric(downwind_propensity_fit$y)  ) ) - 1/nrow(downwind_positive_data)   ) ) -  sum(hatm_0 * ( (1/( 1 - hatpi  )) * (hatpi - 2* as.numeric(downwind_propensity_fit$y) + 1 ) /( sum( (1/(1-hatpi)) * (1 - as.numeric(downwind_propensity_fit$y))  ) )  - 1/nrow(downwind_positive_data)  ) )
  # check3 - check
  # check2 - check
  wrong.sate.aipw = sate.ipw.ma + check2


  # #Checking expression below eq(9) of JRSSA
  # alternative_hatm_0 = lme4::getME(downwind_lmm_fit, 'y') - as.vector(z_mat %*% hatbeta_downwind) * as.numeric(downwind_propensity_fit$y)
  # #alternative_hatm_1 = lme4::getME(downwind_lmm_fit, 'y') + as.vector(z_mat %*% hatbeta_downwind) * (1-as.numeric(downwind_propensity_fit$y))
  # alternative_sate.ipw.ma = mean(alternative_hatm_1) - mean(alternative_hatm_0) + sum(hatw_1 * as.numeric(downwind_propensity_fit$y) * (lme4::getME(downwind_lmm_fit, 'y') - alternative_hatm_1)  ) - sum( hatw_0 * (1 - as.numeric(downwind_propensity_fit$y) ) * (lme4::getME(downwind_lmm_fit, 'y'  ) - alternative_hatm_0)  )
  # # The expression below eq(9) of JRSSA is incorrect
  # alternative_sate.ipw.ma - sate.mb
  # # This is because alternative_sate.ipw.ma is equal to sum_{i=1}^{n} \hat{lambda}_i I_i / n, instead of sum_{i=1}^{n} \hat{lambda}_i I_i / sum_{i=1}^{n} I_i - see its derivation in Ipad's 'Relationship between IPW-MA and MB'
  # alternative_sate.ipw.ma - mean(as.vector(z_mat %*% hatbeta_downwind))


  #Checking eq(10) of JRSSA - can't really check since it involves y_{0i} which is unobserved
  # x_mat = x_z_mat[,which(colnames(x_z_mat) %in%  c('(Intercept)',x_downwind_name) )]
  # hate = lme4::getME(downwind_lmm_fit, 'y') - as.numeric(downwind_propensity_fit$y)*as.vector(z_mat %*% hatbeta_downwind) - as.vector(x_mat %*% hatalpha_downwind)
  # haty_0 = as.vector(x_mat %*% hatalpha_downwind) + hate
  # #This expression is incorrect since we dont really have y_{01}, so we are incorrectly using y_i here
  # R_0 = lme4::getME(downwind_lmm_fit, 'y') - haty_0
  # eq10 = sum(hatw_1 * as.numeric(downwind_propensity_fit$y) * as.vector(z_mat %*% hatbeta_downwind) ) - sum( hatw_1 *  R_0) - sum( hatw_0 * (1 - as.numeric(downwind_propensity_fit$y)) * R_0 )
  # alternative_sate.ipw.l = sum( hatw_1 * as.numeric(downwind_propensity_fit$y) * lme4::getME(downwind_lmm_fit, 'y') ) - sum(hatw_0 * (1-as.numeric(downwind_propensity_fit$y)) * lme4::getME(downwind_lmm_fit, 'y')) - sum(hatw_1 * haty_0) + sum(hatw_0 * (1-as.numeric(downwind_propensity_fit$y)) * haty_0)


  return(list(
    estimates = list(sate.mb = sate.mb,
                     sate.ipw = sate.ipw,
                     sate.ipw.l = sate.ipw.l,
                     sate.ipw.ma = sate.ipw.ma,
                     sate.aipw = sate.aipw,
                     wrong.sate.aipw = wrong.sate.aipw),
    fitted_models = list(
      downwind_propensity_fit = downwind_propensity_fit,
      downwind_positive_target_lmm_fit = downwind_positive_target_lmm_fit,
      downwind_positive_control_lmm_fit = downwind_positive_control_lmm_fit
    )
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

#TODO: Consider to add another variation of 'PREB2' and 'REB2' for adjusting downwind_lmm_fit's fixef as well as random effects, and plug these corrected estimates to compute hatattr and hatsate, rather than directly centering hatattr and hatsate
#TODO: Consider adding a parallelization option
bootstrap_downwind = function(B_bootstrap, bootstrap_type, bootstrap_zero, positive_prob_threshold = NULL, discretize_rain, winsorize_individual_rain, winsorize_total_rain,
                              ori_data, downwind, ori_positive, rain_col_name, downwind_target_subset, downwind_control_subset, ori_fitted_models,
                              downwind_lmm_formula, attr_type, x_downwind_name, target_only,
                              downwind_propensity_formula,
                              ori_attr_est, ori_sate_est){

  if(!bootstrap_type %in% c('REB0','REB1','REB2','PREB0','PREB1','PREB2','MREB1')){
    stop("bootstrap_type should be one of c('REB0','REB1','REB2','PREB0','PREB1','PREB2','MREB1')")
  }



  #Pre-compute hat{P}_{it} which will be used later on to simulate zero vs non-zero rainfall events
  if(bootstrap_zero){
    init_downwind_positive_prob = predict(ori_fitted_models$downwind_logistic_fit,type = 'response')
    if(!is.null(positive_prob_threshold)){
      downwind_positive_prob = init_downwind_positive_prob
      downwind_positive_prob[downwind_positive_prob < positive_prob_threshold] = 0
      downwind_positive_prob[downwind_positive_prob >= positive_prob_threshold] = downwind_positive_prob[downwind_positive_prob >= positive_prob_threshold]* sum(init_downwind_positive_prob)/ sum(downwind_positive_prob[downwind_positive_prob >= positive_prob_threshold])
    }else{
      downwind_positive_prob = init_downwind_positive_prob
    }
  }


  r_vec <- lme4::getME(ori_fitted_models$downwind_lmm_fit, 'y')  - predict(ori_fitted_models$downwind_lmm_fit, re.form = NA)
  group_name = names(lme4::getME(ori_fitted_models$downwind_lmm_fit, "flist"))
  ori_downwind_positive_group = ori_data[downwind & ori_positive , group_name]
  ori_downwind_positive_group_label = unique(ori_downwind_positive_group)
  ori_D_groups =  length(ori_downwind_positive_group_label)

  if(bootstrap_type %in% c('REB0', 'REB2', 'PREB0', 'PREB2')){
    hat.u = sapply(ori_downwind_positive_group_label, function(x){
      mean(r_vec[ori_downwind_positive_group == x])
    })


    hat.e <- rep(0,sum(downwind & ori_positive))
    for(h in 1:ori_D_groups){
      hat.e[ori_downwind_positive_group==ori_downwind_positive_group_label[h]] <- r_vec[ori_downwind_positive_group==ori_downwind_positive_group_label[h]]- hat.u[h]
    }

    final.hat.u = hat.u
    final.hat.e = hat.e
  }

  if(bootstrap_type %in%  c('REB1','PREB1', 'MREB1')){
    hat.u = sapply(ori_downwind_positive_group_label, function(x){
      mean(r_vec[ori_downwind_positive_group == x])
    })

    hat.e <- rep(0,sum(downwind & ori_positive))
    for(h in 1:ori_D_groups){
      hat.e[ori_downwind_positive_group==ori_downwind_positive_group_label[h]] <- r_vec[ori_downwind_positive_group==ori_downwind_positive_group_label[h]]- hat.u[h]
    }


    vc_df <- as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))
    hatsigma2.u = vc_df[vc_df$grp == group_name, "vcov"]
    hat.u.c = hat.u - mean(hat.u)
    if(bootstrap_type == 'REB1'){
      hat.u.cs = ( sqrt(hatsigma2.u) / sqrt( mean(hat.u^2)  )    ) * hat.u.c
    }

    if(bootstrap_type %in% c('PREB1','MREB1')){
      hat.u.cs = ( sqrt(hatsigma2.u) / sqrt( mean(hat.u.c^2)  )    ) * hat.u.c
    }


    hatsigma2.e = vc_df[vc_df$grp == "Residual", "vcov"]
    if(bootstrap_type %in% c('REB1','PREB1')){
      hat.e.s = ( sqrt(hatsigma2.e) / sqrt( mean(hat.e^2)  )    ) * hat.e
    }
    if(bootstrap_type == 'MREB1'){
      ni_vec = sapply(ori_downwind_positive_group, FUN = function(x){sum(ori_downwind_positive_group == x)})
      hat.e.s = ( sqrt(hatsigma2.e) / sqrt( sum( (1/ori_D_groups) * (1/ni_vec) * (hat.e^2)  )  )    ) * hat.e
    }


    final.hat.u = hat.u.cs
    final.hat.e = hat.e.s
  }

  #
  if(bootstrap_type %in% c('REB0','REB1','REB2','MREB1')){
    cluster_sample_prob = rep(1/ori_D_groups, ori_D_groups)
  }

  if(bootstrap_type %in% c('PREB0', 'PREB1', 'PREB2')){
    cluster_sample_prob = sapply(ori_downwind_positive_group_label, function(x){ mean(ori_downwind_positive_group  == x) })
  }



  #Also, we extract ori_downwind_data since our following simulated 'b_downwind_positive_vec' should be of length N_Downwind instead of N
  ori_downwind_data = ori_data[downwind,]
  num_downwind = sum(downwind)

  #To store bootstrapped estimates for attribution and SATE
  bootstrap_attr_matrix = matrix(data = NA, nrow = B_bootstrap, ncol = 2, dimnames = list(NULL, c('apo','apl')))
  bootstrap_sate_matrix = matrix(data = NA, nrow = B_bootstrap, ncol = 5, dimnames = list(NULL, c('sate.mb','sate.ipw','sate.ipw.l','sate.ipw.ma','sate.aipw')))
  bootstrap_downwind_lmm_param_matrix = matrix(data = NA, nrow = B_bootstrap,
                                               ncol = length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))[,'vcov']),
                                               dimnames = list(NULL, c(names(lme4::fixef(ori_fitted_models$downwind_lmm_fit)), paste0('VarComponent_',as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))[,'grp'] ) )))
  if(bootstrap_zero){bootstrap_downwind_logistic_param_matrix = matrix(data = NA, nrow = B_bootstrap,
                                                                       ncol = length(coef(ori_fitted_models$downwind_logistic_fit)),
                                                                       dimnames = list(NULL, names(coef(ori_fitted_models$downwind_logistic_fit))))}
  bootstrap_downwind_propensity_param_matrix = matrix(data = NA, nrow = B_bootstrap,
                                                      ncol = length(coef(ori_fitted_models$downwind_propensity_fit)),
                                                      dimnames = list(NULL, names(coef(ori_fitted_models$downwind_propensity_fit))))
  bootstrap_downwind_positive_target_lmm_param_matrix = matrix(data = NA, nrow = B_bootstrap,
                                                               ncol = length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'vcov']),
                                                               dimnames = list(NULL, c(names(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)), paste0('VarComponent_',as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'grp'] ) )))
  bootstrap_downwind_positive_control_lmm_param_matrix = matrix(data = NA, nrow = B_bootstrap,
                                                                ncol = length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'vcov']),
                                                                dimnames = list(NULL, c(names(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)), paste0('VarComponent_',as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'grp'] ) )))

  bootstrap_downwind_response_matrix = matrix(data = NA, nrow = B_bootstrap, ncol = num_downwind)

  z_downwind_name = setdiff(names(lme4::fixef(ori_fitted_models$downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
  downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)

  #browser()
  for(b in 1:B_bootstrap){
    tryCatch({

      if(bootstrap_zero){
        b_downwind_positive = (runif(num_downwind, min = 0, max = 1) < downwind_positive_prob)
        b_downwind_logistic_fit = glm(b_downwind_positive ~ model.matrix(ori_fitted_models$downwind_logistic_fit$formula, data = ori_data[downwind,]) - 1, data = ori_data[downwind,], family = 'binomial')
        bootstrap_downwind_logistic_param_matrix[b,] = coef(b_downwind_logistic_fit)
      }else{
        b_downwind_positive = ori_positive[downwind]
      }

      b_downwind_positive_data = ori_downwind_data[b_downwind_positive,]
      b_num_downwind_positive = sum(b_downwind_positive)
      b_downwind_positive_group = b_downwind_positive_data[,group_name]
      b_downwind_positive_group_label = unique(b_downwind_positive_group)
      b_D_groups =  length(b_downwind_positive_group_label)

      b_fitted = predict(ori_fitted_models$downwind_lmm_fit, newdata = b_downwind_positive_data, re.form = NA)

      b_y <- rep(NA, b_num_downwind_positive)

      b_donor_group_label <- sample(x = ori_downwind_positive_group_label,
                                    size = b_D_groups,
                                    replace = T,
                                    prob = cluster_sample_prob)

      b_u = sample(x = final.hat.u,  size= b_D_groups, replace=T)



      for(h in 1:b_D_groups){
        target.units = (1:b_num_downwind_positive)[b_downwind_positive_group == b_downwind_positive_group_label[h] ]
        donor.units = (1:sum(downwind & ori_positive))[ori_downwind_positive_group == b_donor_group_label[h]]
        if(length(donor.units) > 1){
          donating.units = sample(x=donor.units, size = length(target.units), replace = T)
        }else{
          donating.units = rep(donor.units, length(target.units))
        }


        b_y[target.units] <- b_fitted[target.units] + b_u[h] + final.hat.e[donating.units]
      }

      #Perform (optional) adjustment of raw rainfall
      b_raw_y = exp(b_y)

      if(discretize_rain){
        b_raw_y[b_raw_y<0.3] <- 0.2
        b_raw_y[(b_raw_y>0.3)&(b_raw_y<0.5)] <- 0.4
        b_raw_y[(b_raw_y>0.5)&(b_raw_y<0.7)] <- 0.6
        b_raw_y[(b_raw_y>0.7)&(b_raw_y<0.9)] <- 0.8
      }

      if(winsorize_individual_rain){
        b_raw_y[b_raw_y>175] <- 100+75*runif(n=sum(b_raw_y>175))
      }

      if(winsorize_total_rain){
        if(sum(b_raw_y)<6000 | sum(b_raw_y)>60000){
          b_raw_y <- b_raw_y*(runif(n=1,min=6000,max=60000))/sum(b_raw_y)
        }
      }

      b_y = log(b_raw_y)

      #This part reconstructs the bootstrapped raw rain, depending on whether there are any offset terms specified in downwind_lmm_formula
      if(length(formula.tools::lhs.vars(downwind_lmm_formula)) == 1){
        b_downwind_positive_data[,rain_col_name] = exp(b_y)
      }

      if(length(formula.tools::lhs.vars(downwind_lmm_formula)) > 1){
        all_offset_terms = formula.tools::lhs.vars(downwind_lmm_formula)[2:length(formula.tools::lhs.vars(downwind_lmm_formula))]
        if(length(all_offset_terms) > 1){
          b_downwind_positive_data[,rain_col_name] = exp(b_y + apply(b_downwind_positive_data[,all_offset_terms],1,sum))
        }

        if(length(all_offset_terms) == 1){
          b_downwind_positive_data[,rain_col_name] = exp(b_y + as.vector(b_downwind_positive_data[,all_offset_terms]))
        }
      }


      #Create a new column 'bootstrapped_y' instead of replacing the LogRain column to accommodate for the case of having LogRain - natural_pred on the RHS of downwind_lmm_formula
      #In this case, we are essentially creating bootstrapped_y = LogRain* - natural_pred, where LogRain* = natural_pred + Xhatbeta + u* + e*
      b_downwind_positive_data$bootstrapped_y = b_y
      b_downwind_lmm_fit = lme4::lmer(update.formula(downwind_lmm_formula, bootstrapped_y ~ . ),
                                      data = b_downwind_positive_data)

      bootstrap_downwind_lmm_param_matrix[b,] = c(lme4::fixef(b_downwind_lmm_fit),
                                                  as.data.frame(lme4::VarCorr(b_downwind_lmm_fit))[,'vcov'])



      bootstrap_downwind_response_matrix[b, b_downwind_positive] = b_y
      bootstrap_downwind_response_matrix[b, !b_downwind_positive] = NA


      # b_downwind_positive_target = b_downwind_positive_data$Gauge.Day.Type == 'Target'
      # b_downwind_positive_control = b_downwind_positive_data$Gauge.Day.Type == 'Control'
      #browser()
      b_downwind_positive_target = rlang::eval_tidy(rlang::enquo(downwind_target_subset), data = b_downwind_positive_data)
      b_downwind_positive_control = rlang::eval_tidy(rlang::enquo(downwind_control_subset), data = b_downwind_positive_data)


      b_hatattr = attr_est(attr_type, b_downwind_positive_data, rain_col_name, b_downwind_positive_target, b_downwind_positive_control,
                           x_downwind_name, target_only, downwind_lmm_fit = b_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

      #Need to also update the downwind_separate_formula here, since the bootstrapped 2nd stage response is now stored in the column 'bootstrapped_y'
      b_hatsate = sate_est(b_downwind_positive_data, b_downwind_positive_target, b_downwind_positive_control, downwind_propensity_formula, update.formula(downwind_separate_formula, bootstrapped_y ~ . ),
                           x_downwind_name, downwind_lmm_fit = b_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

      bootstrap_downwind_propensity_param_matrix[b,] = coef(b_hatsate$fitted_models$downwind_propensity_fit)
      bootstrap_downwind_positive_target_lmm_param_matrix[b,] = c(lme4::fixef(b_hatsate$fitted_models$downwind_positive_target_lmm_fit),
                                                                  as.data.frame(lme4::VarCorr(b_hatsate$fitted_models$downwind_positive_target_lmm_fit))[,'vcov'])
      bootstrap_downwind_positive_control_lmm_param_matrix[b,] = c(lme4::fixef(b_hatsate$fitted_models$downwind_positive_control_lmm_fit),
                                                                   as.data.frame(lme4::VarCorr(b_hatsate$fitted_models$downwind_positive_control_lmm_fit))[,'vcov'])


      #TODO: try to save tor and tlr as well, to understand why the bootstrap plots for apo and apl are always the same shape
      #Maybe can also go back and look at previous plots to see if we always have same shape for apo and apl
      bootstrap_attr_matrix[b,] = c(b_hatattr$apo, b_hatattr$apl)
      bootstrap_sate_matrix[b,] = c(b_hatsate$estimates$sate.mb, b_hatsate$estimates$sate.ipw, b_hatsate$estimates$sate.ipw.l, b_hatsate$estimates$sate.ipw.ma, b_hatsate$estimates$sate.aipw)




    },error=function(e){cat(b,"th","Bootstrap Run Skipped due to ERROR :",conditionMessage(e), "\n")})
  }


  #browser()

  if(bootstrap_type %in% c('REB2','PREB2')){

    bootstrap_downwind_lmm_param_matrix[,(length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + 1): (length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))[,'vcov']) )] =
      adjust_bootstrap_var_components(bootstrapped_var_components = bootstrap_downwind_lmm_param_matrix[,(length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + 1): (length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))[,'vcov']) )])

    bootstrap_downwind_positive_target_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'vcov'])   )] =
      adjust_bootstrap_var_components(bootstrapped_var_components = bootstrap_downwind_positive_target_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'vcov'])   )])

    bootstrap_downwind_positive_control_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'vcov'])   )] =
      adjust_bootstrap_var_components(bootstrapped_var_components = bootstrap_downwind_positive_control_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'vcov'])   )])


    #Perform correction based on the original estimates for each bootstrapped estimates matrix
    bootstrap_attr_matrix = bootstrap_attr_matrix + matrix( data = rep(ori_attr_est - apply(bootstrap_attr_matrix,2, function(x){mean(x, na.rm = T)}), B_bootstrap),
                                                            nrow = B_bootstrap, ncol = ncol(bootstrap_attr_matrix), byrow = TRUE)

    bootstrap_sate_matrix = bootstrap_sate_matrix + matrix( data = rep(ori_sate_est - apply(bootstrap_sate_matrix,2, function(x){mean(x, na.rm = T)}), B_bootstrap),
                                                            nrow = B_bootstrap, ncol = ncol(bootstrap_sate_matrix), byrow = TRUE)

    bootstrap_downwind_lmm_param_matrix[,1: length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) ] =
      bootstrap_downwind_lmm_param_matrix[,1: length(lme4::fixef(ori_fitted_models$downwind_lmm_fit))] +
      matrix( data = rep( lme4::fixef(ori_fitted_models$downwind_lmm_fit) - apply(bootstrap_downwind_lmm_param_matrix[,1: length(lme4::fixef(ori_fitted_models$downwind_lmm_fit))],2, function(x){mean(x, na.rm = T)}), B_bootstrap),
              nrow = B_bootstrap, ncol = ncol(bootstrap_downwind_lmm_param_matrix[,1: length(lme4::fixef(ori_fitted_models$downwind_lmm_fit))]), byrow = TRUE)


    bootstrap_downwind_lmm_param_matrix[,(length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + 1): (length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))[,'vcov']) )] =
      bootstrap_downwind_lmm_param_matrix[,(length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + 1): (length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))[,'vcov']) )] *
      matrix( data = rep( as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))[,'vcov'] / apply(bootstrap_downwind_lmm_param_matrix[,(length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + 1): (length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))[,'vcov']) )], 2, function(x){mean(x,na.rm = T)}), B_bootstrap  ),
              nrow = B_bootstrap, ncol = ncol(bootstrap_downwind_lmm_param_matrix[,(length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + 1): (length(lme4::fixef(ori_fitted_models$downwind_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))[,'vcov']) )]), byrow = TRUE)

    if(bootstrap_zero){
      bootstrap_downwind_logistic_param_matrix = bootstrap_downwind_logistic_param_matrix + matrix(data = rep( coef(ori_fitted_models$downwind_logistic_fit) - apply(bootstrap_downwind_logistic_param_matrix,2, function(x){mean(x, na.rm = T)}), B_bootstrap),
                                                                                                   nrow = B_bootstrap, ncol = ncol(bootstrap_downwind_logistic_param_matrix), byrow = TRUE)
    }

    bootstrap_downwind_propensity_param_matrix = bootstrap_downwind_propensity_param_matrix + matrix(data = rep( coef(ori_fitted_models$downwind_propensity_fit) - apply(bootstrap_downwind_propensity_param_matrix,2, function(x){mean(x, na.rm = T)}), B_bootstrap),
                                                                                                     nrow = B_bootstrap, ncol = ncol(bootstrap_downwind_propensity_param_matrix), byrow = TRUE)

    bootstrap_downwind_positive_target_lmm_param_matrix[, 1: length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit))] =
      bootstrap_downwind_positive_target_lmm_param_matrix[, 1: length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit))] +
      matrix( data = rep( lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit) - apply(bootstrap_downwind_positive_target_lmm_param_matrix[, 1: length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit))], 2, function(x){mean(x, na.rm = T)}), B_bootstrap),
              nrow = B_bootstrap, ncol = ncol(bootstrap_downwind_positive_target_lmm_param_matrix[, 1: length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit))]), byrow = TRUE)

    bootstrap_downwind_positive_target_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'vcov'])   )] =
      bootstrap_downwind_positive_target_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'vcov'])   )] *
      matrix( data = rep( as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'vcov']  / apply(bootstrap_downwind_positive_target_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'vcov'])   )], 2, function(x){mean(x,na.rm=T)}), B_bootstrap),
              nrow = B_bootstrap, ncol = ncol(bootstrap_downwind_positive_target_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'vcov'])   )]), byrow = TRUE)

    bootstrap_downwind_positive_control_lmm_param_matrix[, 1: length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit))] =
      bootstrap_downwind_positive_control_lmm_param_matrix[, 1: length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit))] +
      matrix( data = rep( lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit) - apply(bootstrap_downwind_positive_control_lmm_param_matrix[, 1: length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit))], 2, function(x){mean(x, na.rm = T)}), B_bootstrap),
              nrow = B_bootstrap, ncol = ncol(bootstrap_downwind_positive_control_lmm_param_matrix[, 1: length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit))]), byrow = TRUE)

    bootstrap_downwind_positive_control_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'vcov'])   )] =
      bootstrap_downwind_positive_control_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'vcov'])   )] *
      matrix( data = rep( as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'vcov']  / apply(bootstrap_downwind_positive_control_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'vcov'])   )], 2, function(x){mean(x,na.rm=T)}), B_bootstrap),
              nrow = B_bootstrap, ncol = ncol(bootstrap_downwind_positive_control_lmm_param_matrix[,( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + 1 ): ( length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'vcov'])   )]), byrow = TRUE)

  }

  if(bootstrap_zero){
    return(list(
      hatattr = bootstrap_attr_matrix,
      hatsate = bootstrap_sate_matrix,
      downwind_lmm_param = bootstrap_downwind_lmm_param_matrix,
      downwind_logistic_param = bootstrap_downwind_logistic_param_matrix,
      downwind_propensity_param = bootstrap_downwind_propensity_param_matrix,
      downwind_positive_target_lmm_param = bootstrap_downwind_positive_target_lmm_param_matrix,
      downwind_positive_control_lmm_param = bootstrap_downwind_positive_control_lmm_param_matrix,
      downwind_response = bootstrap_downwind_response_matrix
    ))
  }else{
    return(list(
      hatattr = bootstrap_attr_matrix,
      hatsate = bootstrap_sate_matrix,
      downwind_lmm_param = bootstrap_downwind_lmm_param_matrix,
      downwind_logistic_param = NULL,
      downwind_propensity_param = bootstrap_downwind_propensity_param_matrix,
      downwind_positive_target_lmm_param = bootstrap_downwind_positive_target_lmm_param_matrix,
      downwind_positive_control_lmm_param = bootstrap_downwind_positive_control_lmm_param_matrix,
      downwind_response = bootstrap_downwind_response_matrix
    ))
  }
}

#
adjust_bootstrap_var_components = function(bootstrapped_var_components){
  if(ncol(bootstrapped_var_components)!= 2){
    stop('There should be only 2 columns in bootstrapped_var_components corresponding to the variance components of cluster and residuals')
  }


  L.mat.b <- log(bootstrapped_var_components)
  mu.me <- apply(L.mat.b,2,mean)
  C.mat.b <- cov(L.mat.b)
  su.se <- sqrt(diag(C.mat.b))

  temp <- eigen(solve(C.mat.b),symmetric=T)
  C.mat.b.neg.half <- temp$vectors%*%diag(sqrt(temp$values))
  M.mat.b <- cbind(rep(mu.me[1],nrow(bootstrapped_var_components)),rep(mu.me[2],nrow(bootstrapped_var_components)))

  Sbmod <- (L.mat.b- M.mat.b) %*% C.mat.b.neg.half
  Sbmod[,1] <- Sbmod[,1]* su.se[1]
  Sbmod[,2] <- Sbmod[,2]* su.se[2]

  output <- exp(M.mat.b + Sbmod)
  return(output)
}

#
bootstrap_p_value = function(bootstrap_result){
  if(is.null(bootstrap_result)){
    return(NULL)
  }else{
    return(apply(bootstrap_result,2, function(x){mean(x < 0, na.rm = T)}))
  }
}


bootstrap_CI = function(bootstrap_result, level){
  if(is.null(bootstrap_result)){
    return(NULL)
  }else{
    return(t(apply(bootstrap_result,2, function(x){quantile(x, probs = c( (1-level)/2, 1 - (1-level)/2  ), na.rm = T)})))
  }
}



bootstrap_plot = function(bootstrap_result, ori_est){
  num_var = ncol(bootstrap_result)
  output_plot_list = lapply(1:num_var, function(i){
    ggplot2::ggplot() +
      ggplot2::geom_density(ggplot2::aes(x = bootstrap_result[,i])) +
      ggplot2::geom_vline(xintercept = ori_est[i]) +
      ggplot2::geom_vline(xintercept = 0,linetype="dotted") +
      ggplot2::ggtitle(colnames(bootstrap_result)[i]) +
      ggplot2::xlab('Bootstrapped Values') +
      ggplot2::ylab('Bootstrap Density') +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  })
  names(output_plot_list) = colnames(bootstrap_result)
  return(output_plot_list)
}

#

# perm_HOp_gaugeday = dplyr::left_join(oman[,c('Year','TrialDay')], ionizer_operation, by = c('Year','TrialDay'))
# perm_HOp = perm_HOp_gaugeday[, -which(colnames(perm_HOp_gaugeday) %in% c('Year','TrialDay') )]
#
# #Replace Target.H.XX columns with the elementwise product between HDown and the newly permuted HOp matrix
# mean( (perm_HOp * gaugeday_downwind) == oman[,45:54])
#
#
# #
# B_permutation = 5
# permute_between_ionizer = T
# permute_all_ionizers_between_day = T
# permute_between_gaugeday = T
#
# ionizer_operation = ionizer_operation
# gaugeday_downwind = gaugeday_downwind
# year_ionizer_list =
#   list(
#     '2013' = c('H1','H2'),
#     '2014' = c('H1','H2','H3','H4'),
#     '2015' = c('H1','H2','H3','H4','H5','H6'),
#     '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'),
#     '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10'),
#     '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8', 'H9', 'H10')
#   )
#
# data_target_column_names = names(oman[45:54])
# ionizer_operation_year_column_name = 'Year'
# ionizer_operation_day_column_name = 'TrialDay'
# data = asd3_PREB1$data
# downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay)
# downwind_propensity_formula = Gauge.Day.Type == 'Target' ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure
# rain_col_name = 'Rain.Gauge.Measurement'
# attr_type = 'Proposed'
# x_downwind_name = c('Gauge.Elevation', 'natural_pred')
# target_only = FALSE

#TODO: Consider to add parallelization option
permutation_ionizer = function(B_permutation, permute_between_ionizer, permute_all_ionizers_between_day, permute_between_gaugeday,
                               ionizer_operation, gaugeday_downwind, year_ionizer_list,
                               data_target_column_names, ionizer_operation_year_column_name, ionizer_operation_day_column_name,
                               data, downwind_lmm_formula, downwind_propensity_formula,
                               attr_type, x_downwind_name, target_only,
                               rain_col_name){
  #data_target_column_names are the column names of 'data', which correspond to the target indicators of all ionizers
  #ionizer_operation_year_column_name is the column name of 'ionizer_operation' containing the year of each day
  #ionizer_operation_day_column_name is the column name of 'ionizer_operation' containing the day of each observation, which should be the same column name in 'data'
  #Note that data_target_column_names should ahve the same ordering as colnames(ionizer_operation), as well as gaugeday_downwind
  #The rows of gaugeday_downwind should be ordered as the same as the order of rows of gauge-day observations in 'data'
  if(length(unique(ionizer_operation[,ionizer_operation_day_column_name])) != nrow(ionizer_operation)){
    stop('ionizer_operation has more than one rows associated to the same day')
  }

  perm_attr_matrix = matrix(data = NA, nrow = B_permutation, ncol = 2, dimnames = list(NULL, c('apo','apl')))
  perm_sate_matrix = matrix(data = NA, nrow = B_permutation, ncol = 5, dimnames = list(NULL, c('sate.mb','sate.ipw','sate.ipw.l','sate.ipw.ma','sate.aipw')))

  # ionizer_operation_yearlist = lapply(
  #   names(year_ionizer_list), function(x){
  #     ionizer_operation[ionizer_operation[,ionizer_operation_year_column_name] == x, -which(colnames(ionizer_operation) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name) )  ]
  #   }
  # )
  # names(ionizer_operation_yearlist) = names(year_ionizer_list)

  for(b in 1:B_permutation){
    tryCatch({
    perm_data = data



    perm_ionizer_operation_day = ionizer_operation
    if(permute_between_ionizer){
      for(unique_year in unique(perm_ionizer_operation_day[, ionizer_operation_year_column_name])){
        temp = perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))]
        deployed_ionizers = year_ionizer_list[[as.character(unique_year)]]
        perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))] =
          t(apply(temp,1, FUN = function(x){
            x[colnames(temp) %in% deployed_ionizers] = sample(x[colnames(temp) %in% deployed_ionizers])
            return(x)
          }))
      }
    }

    if(permute_all_ionizers_between_day){
      for(unique_year in unique(perm_ionizer_operation_day[, ionizer_operation_year_column_name])){
        temp = perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))]
        perm_ionizer_operation_day[perm_ionizer_operation_day[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_day_column_name,ionizer_operation_year_column_name))] =
          temp[sample(1:nrow(temp)),]
      }
    }


    # perm_ionizer_operation_yearlist = ionizer_operation_yearlist
    # if(permute_between_ionizer){
    #   for(i in 1:length(perm_ionizer_operation_yearlist)){
    #     deployed_ionizers = year_ionizer_list[[names(perm_ionizer_operation_yearlist)[i]]]
    #     perm_ionizer_operation_yearlist[[i]] = t(apply(perm_ionizer_operation_yearlist[[i]], 1, function(x){
    #       x[colnames(x) %in% deployed_ionizers] = sample(x[colnames(x) %in% deployed_ionizers])
    #       return(x)
    #     }))
    #   }
    # }
    #
    # if(permute_all_ionizers_between_day){
    #   for(i in 1:length(perm_ionizer_operation_yearlist)){
    #     perm_ionizer_operation_yearlist[[i]] = perm_ionizer_operation_yearlist[[i]][sample(1:nrow(perm_ionizer_operation_yearlist[[i]])),]
    #   }
    # }
    #
    # perm_ionizer_operation_day = ionizer_operation
    # for(i in 1:length(year_ionizer_list)){
    #   perm_ionizer_operation_day[perm_ionizer_operation_day[,ionizer_operation_year_column_name] == names(year_ionizer_list)[i], -which(colnames(perm_ionizer_operation_day) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name) )  ] = perm_ionizer_operation_yearlist[[i]]
    # }

    perm_ionizer_operation_gaugeday = dplyr::left_join(perm_data[, ionizer_operation_day_column_name, drop = FALSE], perm_ionizer_operation_day, by = ionizer_operation_day_column_name)

    if(permute_between_gaugeday){
      for(unique_year in unique(perm_ionizer_operation_gaugeday[, ionizer_operation_year_column_name])){
        temp = perm_ionizer_operation_gaugeday[perm_ionizer_operation_gaugeday[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_gaugeday) %in% c(ionizer_operation_day_column_name, ionizer_operation_year_column_name))  ]
        perm_ionizer_operation_gaugeday[perm_ionizer_operation_gaugeday[, ionizer_operation_year_column_name] == unique_year, -which(colnames(perm_ionizer_operation_gaugeday) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name))  ] =
          temp[sample(1: nrow(temp)),]
      }
    }

    perm_data[,data_target_column_names] = (perm_ionizer_operation_gaugeday[, -which(colnames(perm_ionizer_operation_gaugeday) %in% c(ionizer_operation_year_column_name, ionizer_operation_day_column_name))] * gaugeday_downwind)

    #Refit downwind LMM to permuted data, and then compute hatattr and hatsate, noting we need to becareful with the definition of downwind_positive_target, downwind_positive_control, and downwind_propensity_
    downwind = apply(gaugeday_downwind,1,sum) > 0
    positive = ( data[,rain_col_name] > 0)
    perm_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = perm_data[downwind & positive,])

    perm_downwind_positive_data = perm_data[downwind & positive,]
    target_vec = apply(perm_data[,data_target_column_names],1, sum) > 0
    nontarget_vec = apply(perm_data[,data_target_column_names],1, sum) ==  0
    perm_downwind_positive_target = target_vec[downwind & positive]
    perm_downwind_positive_control = nontarget_vec[downwind & positive]

    perm_hatattr = attr_est(attr_type, perm_downwind_positive_data, rain_col_name, perm_downwind_positive_target, perm_downwind_positive_control,
                            x_downwind_name, target_only = target_only, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)


    perm_downwind_positive_data$permuted_target_indicator = as.logical(perm_downwind_positive_target)

    if(b == 1){
      z_downwind_name = setdiff(names(lme4::fixef(perm_downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
      perm_downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)
      perm_downwind_propensity_formula = update.formula(downwind_propensity_formula, permuted_target_indicator ~ . )
    }

    perm_hatsate = sate_est(perm_downwind_positive_data, perm_downwind_positive_target, perm_downwind_positive_control, perm_downwind_propensity_formula, perm_downwind_separate_formula,
                            x_downwind_name, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

    perm_attr_matrix[b,] = c(perm_hatattr$apo, perm_hatattr$apl)
    perm_sate_matrix[b,] = c(perm_hatsate$estimates$sate.mb, perm_hatsate$estimates$sate.ipw, perm_hatsate$estimates$sate.ipw.l, perm_hatsate$estimates$sate.ipw.ma, perm_hatsate$estimates$sate.aipw)
    },error=function(e){cat(b,"th","Permutation Run Skipped due to ERROR :",conditionMessage(e), "\n")})
  }


  return(list(
    hatattr = perm_attr_matrix,
    hatsate = perm_sate_matrix
  ))
}

permutation_p_value = function(permutation_result, ori_est){
  return(sapply(1:ncol(permutation_result), FUN = function(i){
    mean(permutation_result[,i] >= ori_est[i], na.rm = T )
  }))
}

permutation_plot = function(permutation_result, ori_est){
  num_var = ncol(permutation_result)
  output_plot_list = lapply(1:num_var, function(i){
    ggplot2::ggplot() +
      ggplot2::geom_density(ggplot2::aes(x = permutation_result[,i])) +
      ggplot2::geom_vline(xintercept = ori_est[i]) +
      ggplot2::ggtitle(colnames(permutation_result)[i]) +
      ggplot2::xlab('Permuted Values') +
      ggplot2::ylab('Permutation Density') +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  })
  names(output_plot_list) = colnames(permutation_result)
  return(output_plot_list)
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
#                  downwind_control_subset = Gauge.Day.Type == 'Control',
#                  attr_type = 'Proposed',
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
#                              downwind_control_subset = Gauge.Day.Type == 'Control',
#                              attr_type = 'Proposed',
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
# attr_type = 'Proposed'

# asd2  = rain_attr(data = oman,
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
#                  downwind_control_subset = Gauge.Day.Type == 'Control',
#                  attr_type = 'Proposed',
#                  x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                  target_only = FALSE,
#                  bootstrap =T,
#                  bootstrap_option = list(B_bootstrap = 3,
#                                          bootstrap_type = 'REB1',
#                                          bootstrap_zero = T,
#                                          positive_prob_threshold = NULL,
#                                          discretize_rain = T,
#                                          winsorize_individual_rain = T,
#                                          winsorize_total_rain = T,
#                                          CI_level = 0.95
#                                          )
#                  )
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
# asd2$bootstrap_result$downwind_response[1:3,1:10]
# apply(asd2$bootstrap_result$downwind_response, 1 , function(x){sum(!is.na(x))})
#
# #Trying to replicate PREB-1 bootstrap paper results:
# set.seed(123)
# asd3_PREB1  = rain_attr(data = oman,
#                   upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                   instr_pred_name = 'natural_pred',
#                   instr_pred_type = 'Unconditional',
#                   downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                   downwind_logistic_formula = NULL,
#                   downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                   rain_col_name = 'Rain.Gauge.Measurement',
#                   upwind_subset = Gauge.Day.Type == 'Upwind',
#                   downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                   downwind_target_subset = Gauge.Day.Type == 'Target',
#                   downwind_control_subset = Gauge.Day.Type == 'Control',
#                   attr_type = 'No',
#                   x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                   target_only = FALSE,
#                   bootstrap =T,
#                   bootstrap_option = list(B_bootstrap = 3,
#                                           bootstrap_type = 'PREB1',
#                                           bootstrap_zero = F,
#                                           positive_prob_threshold = NULL,
#                                           discretize_rain = F,
#                                           winsorize_individual_rain = F,
#                                           winsorize_total_rain = F,
#                                           CI_level = 0.95
#                   )
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
#                         downwind_control_subset = Gauge.Day.Type == 'Control',
#                         attr_type = 'No',
#                         x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                         target_only = FALSE,
#                         bootstrap =T,
#                         bootstrap_option = list(B_bootstrap = 10,
#                                                 bootstrap_type = 'PREB2',
#                                                 bootstrap_zero = F,
#                                                 positive_prob_threshold = NULL,
#                                                 discretize_rain = F,
#                                                 winsorize_individual_rain = F,
#                                                 winsorize_total_rain = F,
#                                                 CI_level = 0.95
#                         )
# )
# apply(asd3_PREB2$bootstrap_result$downwind_positive_target_lmm_param,2,mean)
# asd3_PREB2$all_fitted_models$downwind_positive_target_lmm_fit
#
#
# #Testing REB2
# set.seed(123)
# asd3_REB2  = rain_attr(data = oman,
#                         upwind_lmm_formula = LogRain ~  Gauge.Elevation + Steering.Wind.Speed + Total.Totals + PC2.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure + (1|TrialDay),
#                         instr_pred_name = 'natural_pred',
#                         instr_pred_type = 'Unconditional',
#                         downwind_lmm_formula = LogRain ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02 + (1|TrialDay),
#                         downwind_logistic_formula = (Rain.Gauge.Measurement > 0) ~ Gauge.Elevation + natural_pred  + Target.H.01 + Target.H.02 + Target.H.03 + Target.H.04 + Target.H.05 + Target.H.06 + Target.H.07 + Target.H.08 + Target.H.09 + Target.H.10 + Gauge.Elevation:Target.H.01 + Gauge.Elevation:Target.H.02,
#                         downwind_propensity_formula = (Gauge.Day.Type == 'Target') ~ Total.Totals + PC1.Dry.Temperature + PC1.Relative.Humidity + PC1.Ground.Level.Pressure,
#                         rain_col_name = 'Rain.Gauge.Measurement',
#                         upwind_subset = Gauge.Day.Type == 'Upwind',
#                         downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),
#                         downwind_target_subset = Gauge.Day.Type == 'Target',
#                         downwind_control_subset = Gauge.Day.Type == 'Control',
#                         attr_type = 'No',
#                         x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                         target_only = FALSE,
#                         bootstrap =T,
#                         bootstrap_option = list(B_bootstrap = 10,
#                                                 bootstrap_type = 'REB2',
#                                                 bootstrap_zero = F,
#                                                 positive_prob_threshold = NULL,
#                                                 discretize_rain = F,
#                                                 winsorize_individual_rain = F,
#                                                 winsorize_total_rain = F,
#                                                 CI_level = 0.95
#                         )
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
#                                       downwind_control_subset = Gauge.Day.Type == 'Control',
#                                       attr_type = 'No',
#                                       x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                       target_only = FALSE,
#                                       bootstrap =T,
#                                       bootstrap_option = list(B_bootstrap = 10,
#                                                               bootstrap_type = 'REB2',
#                                                               bootstrap_zero = T,
#                                                               positive_prob_threshold = NULL,
#                                                               discretize_rain = F,
#                                                               winsorize_individual_rain = F,
#                                                               winsorize_total_rain = F,
#                                                               CI_level = 0.95
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
#                                                downwind_control_subset = Gauge.Day.Type == 'Control',
#                                                attr_type = 'No',
#                                                x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                                target_only = FALSE,
#                                                bootstrap =T,
#                                                bootstrap_option = list(B_bootstrap = 10,
#                                                                        bootstrap_type = 'REB2',
#                                                                        bootstrap_zero = T,
#                                                                        positive_prob_threshold = NULL,
#                                                                        discretize_rain = T,
#                                                                        winsorize_individual_rain = T,
#                                                                        winsorize_total_rain = T,
#                                                                        CI_level = 0.95
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
#                     downwind_control_subset = Gauge.Day.Type == 'Control',
#                     attr_type = 'No',
#                     x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                     target_only = FALSE,
#                     bootstrap =F,
#                     bootstrap_option = NULL,
#                     permutation = T,
#                     permutation_option = list(
#                       B_permutation = 5,
#                       permute_between_ionizer = T,
#                       permute_all_ionizers_between_day = T,
#                       permute_between_gaugeday = T,
#                       ionizer_operation = ionizer_operation,
#                       gaugeday_downwind = gaugeday_downwind,
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
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
set.seed(123)
my_perm_result_TT_RayWinsorize = rain_attr(data = oman,
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
                                           downwind_control_subset = Gauge.Day.Type == 'Control',
                                           attr_type = 'Ray Winsorize',
                                           x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
                                           target_only = FALSE,
                                           bootstrap =F,
                                           bootstrap_option = NULL,
                                           permutation = T,
                                           permutation_option = list(
                                             B_permutation = 6,
                                             permute_between_ionizer = T,
                                             permute_all_ionizers_between_day = T,
                                             permute_between_gaugeday = F,
                                             ionizer_operation = ionizer_operation,
                                             gaugeday_downwind = gaugeday_downwind,
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
                                             ionizer_operation_day_column_name = 'TrialDay'
                                           )
)

load('D:/Postdoc/Simulation/Replicate ISR Results/Rdata/permutation_result_Oman_Trial_Data_perm_row_between_gauge_day_F.Rdata')
max(abs(perm_result_TT$perm_attribution_Ray_winsorize_matrix[1:6,c('apo','apl')] - my_perm_result_TT_RayWinsorize$permutation_result$hatattr))


RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
set.seed(123)
my_perm_result_TT_Proposed = rain_attr(data = oman,
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
                                           downwind_control_subset = Gauge.Day.Type == 'Control',
                                           attr_type = 'Proposed',
                                           x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
                                           target_only = FALSE,
                                           bootstrap =F,
                                           bootstrap_option = NULL,
                                           permutation = T,
                                           permutation_option = list(
                                             B_permutation = 6,
                                             permute_between_ionizer = T,
                                             permute_all_ionizers_between_day = T,
                                             permute_between_gaugeday = F,
                                             ionizer_operation = ionizer_operation,
                                             gaugeday_downwind = gaugeday_downwind,
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
                                             ionizer_operation_day_column_name = 'TrialDay'
                                           )
)
max(abs(perm_result_TT$perm_attribution_proposed_matrix[1:6,c('apo','apl')] - my_perm_result_TT_Proposed$permutation_result$hatattr))


#Replicate previosu analysis with permute_between_gaugeday = T
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
set.seed(123)
my_perm_result_TT_RayWinsorize = rain_attr(data = oman,
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
                                           downwind_control_subset = Gauge.Day.Type == 'Control',
                                           attr_type = 'Ray Winsorize',
                                           x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
                                           target_only = FALSE,
                                           bootstrap =F,
                                           bootstrap_option = NULL,
                                           permutation = T,
                                           permutation_option = list(
                                             B_permutation = 6,
                                             permute_between_ionizer = T,
                                             permute_all_ionizers_between_day = T,
                                             permute_between_gaugeday = T,
                                             ionizer_operation = ionizer_operation,
                                             gaugeday_downwind = gaugeday_downwind,
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
                                             ionizer_operation_day_column_name = 'TrialDay'
                                           )
)

load('D:/Postdoc/Simulation/Replicate ISR Results/Rdata/permutation_result_Oman_Trial_Data_perm_row_between_gauge_day_T.Rdata')
max(abs(perm_result_TT$perm_attribution_Ray_winsorize_matrix[1:6,c('apo','apl')] - my_perm_result_TT_RayWinsorize$permutation_result$hatattr))

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
set.seed(123)
my_perm_result_TT_Proposed = rain_attr(data = oman,
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
                                       downwind_control_subset = Gauge.Day.Type == 'Control',
                                       attr_type = 'Proposed',
                                       x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
                                       target_only = FALSE,
                                       bootstrap =F,
                                       bootstrap_option = NULL,
                                       permutation = T,
                                       permutation_option = list(
                                         B_permutation = 6,
                                         permute_between_ionizer = T,
                                         permute_all_ionizers_between_day = T,
                                         permute_between_gaugeday = T,
                                         ionizer_operation = ionizer_operation,
                                         gaugeday_downwind = gaugeday_downwind,
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
                                         ionizer_operation_day_column_name = 'TrialDay'
                                       )
)
max(abs(perm_result_TT$perm_attribution_proposed_matrix[1:6,c('apo','apl')] - my_perm_result_TT_Proposed$permutation_result$hatattr))

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

#
# set.seed(123)
# sample(1:2)
#
# set.seed(123)
# sample(c('a','b'))

#TODO: Do a full-scale replication of previous permutation analysis, bootstrap analysis in Bootstrap paper, and previous bootstrap analysis (involving bootstrap_zero) to verify the correctness of all functions
