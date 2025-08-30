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
# attr_type = 'Chambers_Chandra'
#
# bootstrap_zero = TRUE



# #data is the full data
# #To use instr_pred or natural_pred as an offset, just specify the response to be LogRain - natural_pred in downwind_lmm_formula
# #To fit first-stage model to Control^#, jut specify upwind_subset = Gauge.Day.Type == 'Control'
# rm(list=ls())
# load('data/oman.rda')
# load('data/ionizer_operation.rda')
# load('data/gaugeday_downwind.rda')



#' @title
#' Attribution and Sample Average Treatment Effect for Rainfall Enhancement Trial Data
#'
#' @description
#' Perform estimation and inference of attribution and sample average treatment effect for rainfall enhancement trial data, based on the two-stage linear mixed model (LMM) approach employed in Chambers et al. (2022).
#'
#' @details
#' This function implements a two-stage modelling procedure via the following steps:
#' \enumerate{
#'  \item Fit an upwind (first stage) LMM using \code{lme4::lmer(upwind_lmm_formula)} to the subset of observations from \code{data} satisfying \code{upwind_subset & positive_subset}. This fitted LMM is used to obtain fitted values (named as \code{instr_pred_name}) to be used in the downwind (second stage) LMM.
#'
#'  \item Fit a downwind (second stage) LMM using  \code{lme4::lmer(downwind_lmm_formula)} to the subset of observations from \code{data} satisfying \code{downwind_subset & positive_subset}.
#'  \deqn{y_{ij} =  x_{ij}^\top \alpha + z_{ij}^\top \beta + u_i + e_{ij}  }
#'  where \eqn{i} indexes day (group) and \eqn{j} indexes gauge (unit within group), \eqn{x_{ij}} is a vector of covariates (with names supplied in \code{x_downwind_name}, including intercept) that are not related to the ionizers (treatment), \eqn{z_{ij}} is a vector of ionizer (treatment) related covariates, \eqn{u_i} are random intercepts, and \eqn{e_{ij}} are error terms.
#' }
#' The fitted values obtained from the upwind (first stage) LMM can either be:
#' \itemize{
#'   \item Included as a covariate on the right-hand side of \code{downwind_lmm_formula}, e.g., \code{instr_pred_name = "natural_pred"} and \code{downwind_lmm_formula = LogRain ~ natural_pred + ...}, or
#'   \item Included as an offset term by subtracting it from the response on the left-hand side, e.g., \code{instr_pred_name = "natural_pred"} and \code{downwind_lmm_formula = LogRain - natural_pred ~ ...}.
#' }
#'
#' \strong{Attribution} \cr
#' Two attribution estimates, namely \code{apo} and \code{apl} are computed based on the estimated fixed effect coefficients \eqn{\hat{\alpha}}, \eqn{\hat{\beta}} and EBLUPs \eqn{\hat{u}_i} from the fitted downwind (second stage) LMM. \code{apo} represents the total increase or decrease in downwind rainfall attributed to the ionizer (treatment) as a proportion of the total amount of observed downwind rainfall., while \code{apl} represents the total increase or decrease in downwind rainfall attributed to the ionizer (treatment) as a proportion of the total expected amount of downwind rainfall without the effect of ionizer (treatment).
#' This function allows for three different ways of estimating \code{apo} and \code{apl} as specified by the argument \code{attr_type}:
#' \describe{
#' \item{\code{Chambers_Chandra}}{Attribution is estimated based on the approach of Chambers et al. (2022), to adjust for back-transformation bias due to the modelling of log-transformed rainfall:
#'     \deqn{
#'     \code{apo} = \sum_{(i,j)} Rain_{ij} [ 1 - \max\{\lambda^{-1} \exp(-z_{ij}^\top \hat{\beta}), 0.5\} ] /  \sum_{(i,j)}Rain_{ij}, \quad
#'     \code{apl} = \sum_{(i,j)} Rain_{ij} [ 1 - \max\{\lambda^{-1} \exp(-z_{ij}^\top \hat{\beta}), 0.5\} ] /  \sum_{(i,j)}Rain_{ij} \max\{\lambda^{-1} \exp(-z_{ij}^\top \hat{\beta}), 0.5\},
#'     }
#'     where the summation is either across all observations satisfying \code{downwind_subset & positive_subset} (when \code{target_only = FALSE}), or across all observations satisfying \code{downwind_target_subset & positive_subset} (when \code{target_only = TRUE}), \eqn{Rain_{ij}} is the observed raw rainfall (contained in the column specified by \code{rain_col_name}),
#'     \deqn{
#'     \lambda = 1 + \frac{\sqrt{ (1+m)^2 + 4(\mu - 1)m  } - (1+m)}{2m}, m = \frac{\hat{V}( x_{ij}^\top \hat{\alpha}  + \hat{u}_i ) }{\hat{V}(z_{ij}^\top \hat{\beta})},
#'     }
#'     with \eqn{\hat{V}(\cdot)} denoting the empirical variance either across all observations satisfying \code{downwind_subset & positive_subset} (when \code{target_only = FALSE}), or across all observations satisfying \code{downwind_target_subset & positive_subset} (when \code{target_only = TRUE}), and
#'     \deqn{
#'     \mu = \frac{1}{N} \sum_{(i,j)} \frac{Rain_{ij}}{\exp( x_{ij}^\top \hat{\alpha} + z_{ij}^\top \hat{\beta} + \hat{u}_i )},
#'     }
#'     and \eqn{N} is either the total number of observations satisfying \code{downwind_subset & positive_subset} (when \code{target_only = FALSE}), or the total number of observations satisfying \code{downwind_target_subset & positive_subset} (when \code{target_only = TRUE}).
#'     When an offset term is included on the LHS of \code{downwind_lmm_formula}, the expressions of \eqn{m} and \eqn{\mu} become
#'     \deqn{
#'     m = \frac{\hat{V}( offset_{ij} + x_{ij}^\top \hat{\alpha}  + \hat{u}_i ) }{\hat{V}(z_{ij}^\top \hat{\beta})}, \mu = \frac{1}{N} \sum_{(i,j)} \frac{Rain_{ij}}{\exp( offset_{ij} + x_{ij}^\top \hat{\alpha} + z_{ij}^\top \hat{\beta} + \hat{u}_i )}.
#'     }
#'
#'   }
#'
#' \item{\code{ThoEtAl}}{Attribution is estimated based on an alternative adjustment using the estimated covariance matrix \eqn{\hat{\Sigma}} of \eqn{\hat{\beta}}.
#'     \deqn{
#'     \code{apo} = \sum_{(i,j)} Rain_{ij} \{ 1 - \exp(z_{ij}^\top \hat{\beta} - 0.5 z_{ij}^\top \hat{\Sigma} z_{ij} ) \}/ \sum_{(i,j)}Rain_{ij}, \quad
#'     \code{apl} = \sum_{(i,j)} Rain_{ij} \{ 1 - \exp(z_{ij}^\top \hat{\beta} - 0.5 z_{ij}^\top \hat{\Sigma} z_{ij} ) \} /  \sum_{(i,j)}Rain_{ij} \exp(-z_{ij}^\top \hat{\beta} - 0.5 z_{ij}^\top \hat{\Sigma} z_{ij} ).
#'     }
#' }
#'
#' \item{\code{No}}{Attribution is estimated based on no adjustment.
#'     \deqn{
#'     \code{apo} = \sum_{(i,j)} Rain_{ij} \{ 1 - \exp(z_{ij}^\top \hat{\beta} ) \}/ \sum_{(i,j)}Rain_{ij}, \quad
#'     \code{apl} = \sum_{(i,j)} Rain_{ij} \{ 1 - \exp(z_{ij}^\top \hat{\beta} ) \} /  \sum_{(i,j)}Rain_{ij} \exp(-z_{ij}^\top \hat{\beta} ).
#'     }
#' }
#' All attribution estimates above implicitly assume that the upwind (first stage) and downwind (second stage) LMMs are modelling the log-transformed rainfall instead of the raw rainfall.
#' These attribution estimates should therefore only be interpreted and used when the LHS of both \code{upwind_lmm_formula} and \code{downwind_lmm_formula} contains the log-transformed rainfall, not the raw rainfall.
#'}
#'
#' \strong{SATE} \cr
#' The computation of SATE estimates involves fitting a downwind (second stage) propensity score model using \code{glm(downwind_propensity_formula, family = "binomial")} to the subset of observations from \code{data} satisfying \code{downwind_subset & positive_subset}, with the response being an indicator \eqn{I_{ij}} for whether each observation is exposed to the ionizer (treatment), i.e., \eqn{I_{ij} = 1} if it satisfies \code{downwind_target_subset & positive_subset}, and \eqn{I_{ij} = 0} if it satisfies \code{downwind_control_subset & positive_subset}.
#' The estimated propensity scores (i.e., fitted values) from this fitted propensity score model, denoted as \eqn{\hat{\pi}_{ij}}, are then used to compute the inverse propensity weights (IPW) \eqn{\hat{w}_{ij,1} = \hat{\pi}_{ij}^{-1} / \sum_{(k,l)}  (\hat{\pi}_{kl}^{-1} I_{kl}) } and \eqn{ \hat{w}_{ij,0} = (1 - \hat{\pi}_{ij})^{-1} / \sum_{(k,l)} \{ (1- \hat{\pi}_{kl})^{-1} (1- I_{kl}) \} }, where the summation is over all observations from \code{data} satisfying \code{downwind_subset & positive_subset}.
#' These IPW weights are then used, together with the estimation results of the downwind (second stage) LMM, to obtain the following five types of SATE estimates:
#' \itemize{
#' \item{\code{sate.mb} \eqn{ = \sum_{(i,j)} I_{ij} (z_{ij}^\top \hat{\beta}) /  \sum_{(i,j)} I_{ij} }.
#' }
#'
#' \item{\code{sate.ipw} \eqn{ = \{ \sum_{(i,j)} \hat{w}_{ij,1} I_{ij} y_{ij} \} - \{ \sum_{(i,j)} \hat{w}_{ij,0} (1-I_{ij}) y_{ij} \}  }, where \eqn{y_{ij}} denote the response variable of the downwind (second stage) LMM that might contain offset term.
#' }
#'
#' \item{\code{sate.ipw.l} \eqn{ = \sum_{(i,j)} \hat{w}_{ij,1} I_{ij} z_{ij}^\top \hat{\beta} }.
#' }
#'
#' \item{\code{sate.ipw.ma} \eqn{ = (N^{-1} \sum_{(i,j)} \hat{m}_{ij,1} ) - (N^{-1} \sum_{(i,j)} \hat{m}_{ij,0} ) + \{ \sum_{(i,j)} \hat{w}_{ij,1} I_{ij} (y_{ij} - \hat{m}_{ij,1}) \} - \{ \sum_{(i,j)} \hat{w}_{ij,0} (1- I_{ij}) (y_{ij} - \hat{m}_{ij,0}) \} }, where \eqn{N} is the total number of observations satisfying \code{downwind_subset & positive_subset}.
#' \eqn{\hat{m}_{ij,1}} are fitted values (using fixed effect only) obtained from the downwind (second stage) treatment-only LMM, based on a modified version of \code{downwind_lmm_formula} whose RHS only contains non-ionizer (non-treatment) related covariates \eqn{x_{ij}} and excludes ionizer (treatment) related covariates \eqn{z_{ij}} fitted to the subset of observations from \code{data} satisfying \code{downwind_target_subset & positive_subset}.
#' Similarly, \eqn{\hat{m}_{ij,0}} are fitted values (using fixed effect only) obtained from the downwind (second stage) control-only LMM, based on the same modified version of \code{downwind_lmm_formula} fitted to another subset of observations from \code{data} satisfying \code{downwind_control_subset & positive_subset}.
#' }
#'
#' \item{\code{sate.aipw} \eqn{ = [\sum_{(i,j)} \hat{w}_{ij,1} \{ I_{ij} y_{ij} - (I_{ij} - \hat{\pi}_{ij}) \hat{m}_{ij,1} \} ] - [\sum_{(i,j)} \hat{w}_{ij,0} \{ (1 - I_{ij}) y_{ij} - (I_{ij} - \hat{\pi}_{ij}) \hat{m}_{ij,0} \} ]  }.
#' }
#' }
#'
#' \strong{Bootstrap and Permutation Inference} \cr
#' This function can also be used to perform bootstrap inference on the attribution and SATE, by setting \code{bootstrap = TRUE} and supplying the relevant bootstrap options using \code{\link{bootstrap_option}()}.
#' For full details of the bootstrap procedure, please see \code{\link{bootstrap_downwind}}. Briefly, the bootstrap is carried out in two levels by conditioning on the upwind (first stage) LMM and its fitted values:
#' \itemize{
#' \item{First level is an optional level that is only carried out when \code{bootstrap_option(bootstrap_zero = TRUE)}. This level considers generating bootstrap samples of rainfall event indicator for the subset of observations satisfying \code{downwind_subset} using the predicted probabilities from the downwind logistic model. This model is fitted using \code{glm(downwind_logistic_formula, family = "binomial")} to the subset of observations from \code{data} satisfying \code{downwind_subset}, with the response being an indicator for whether the observed rainfall is greater than zero, i.e., the indicator is defined by the logical expression \code{positive_subset}.
#'   }
#'
#' \item{Second level generates bootstrap samples of positive rainfall for the subset of observations not only satisfying \code{downwind_subset} but also with the first-level bootstrapped rainfall event indicator being equal to one. When \code{bootstrap_option(bootstrap_zero = FALSE)}, then this level generates bootstrap samples of positive rainfall for the subset of observations satisfying \code{downwind_subset & positive_subset}.
#'   This is done using one of the semiparametric bootstrap methods of Chambers & Chandra (2013) and Tho et al. (2025), which involves the use of marginal residuals from the fitted downwind (second stage) LMM.   }
#' }
#' The above attribution and SATE estimates are then computed based on each bootstrap sample of the positive rainfall, forming their respective bootstrap distributions. This function also provides bootstrap distributions of parameters associated with the downwind LMM (\code{downwind_lmm_formula}), downwind logistic model (\code{downwind_logistic_formula}), downwind propensity score model (\code{downwind_propensity_formula}), downwind treatment-only LMM, and downwind control-only LMM.
#' These bootstrap distributions are then used to compute bootstrap p-values (proportion of bootstrapped estimates that are negative), form bootstrap percentile confidence intervals (with confidence level specified in \code{bootstrap_option$CI_level}), and generate their respective plots.
#' It is worth noting that the entire bootstrap procedure (including rainfall resampling, model fitting, and parameter estimation) can be run in parallel by setting \code{bootstrap_option$bootstrap_parallel = TRUE}, using \code{bootstrap_option$bootstrap_parallel_num_worker} workers.
#'
#' Finally, this function enables permutation-based inference on the attribution and SATE, by setting \code{permutation = TRUE} and supplying the relevant permutation options using \code{\link{permutation_option}()}.
#' For full details of the permutation-based procedure, please see \code{\link{permutation_ionizer}}. In short, the permutation-based procedure involves randomly permuting the operating schedules of the ionizers (treatment) and re-estimating the attribution and SATE based on the permuted data, from which permutations distributions of attribution and SATE estimates are formed.
#' These permutation distributions are used to compute permutation p-values (proportion of permuted estimates that are greater than the observed estimates) and generate their respective plots.
#'
#'
#' @param data A data frame containing the variables named in \code{upwind_lmm_formula}, \code{downwind_lmm_formula}, \code{downwind_logistic_formula} (if specified), and \code{downwind_propensity_formula}.
#' It should also contain variables named in \code{rain_col_name}, \code{upwind_subset}, \code{downwind_subset}, \code{downwind_target_subset}, and \code{downwind_control_subset}.
#' @param upwind_lmm_formula A two sided linear formula object to be used in \link[lme4]{lmer}, describing both the fixed-effects and random intercept part of the upwind (first stage) LMM.
#' @param instr_pred_name A character string to store the variable name of the fitted values generated from the upwind (first stage) LMM.
#' @param instr_pred_type Type of fitted values generated from the upwind (first stage) LMM.. If "Unconditional" the fitted values equal to only the estimated fixed effects. If "Conditional" the fitted values equal to the sum of estimated fixed effects and EBLUPs of random intercepts.
#' @param downwind_lmm_formula A two sided linear formula object to be used in \link[lme4]{lmer}, describing both the fixed-effects and random intercept part of the downwind (second stage) LMM. This formula should contain the variable name specified in \code{instr_pred_name}.
#' @param downwind_logistic_formula An optional two sided linear formula object to be used in \code{\link{stats}{glm}} with \code{family = "binomial"}, for fitting a logistic model to the indicators of rainfall event. This only needs to be specified when \code{bootstrap = TRUE} and \code{bootstrap_option$bootstrap_zero = TRUE}.
#' @param downwind_propensity_formula A two sided linear formula object to be used in \code{\link{stats}{glm}} with \code{family = "binomial"}, for fitting a propensity score model to the treatment indicators of downwind (second stage) observations.
#' @param rain_col_name A character string that refers to the column name of the raw scale rainfall in \code{data}.
#' @param upwind_subset A logical expression used to extract the relevant subset of observations from \code{data} to be used in the upwind (first stage) LMM fitting. For example, \code{Gauge.Day.Type == "Upwind"}.
#' @param downwind_subset A logical expression used to extract the relevant subset of observations from \code{data} to be used in the downwind (second stage) LMM fitting. For example, \code{Gauge.Day.Type \%in\% c("Target","Control")}.
#' @param downwind_target_subset A logical expression used to extract the relevant subset of downwind (second stage) observations from \code{data} that were exposed to treatment (operating ionizers). For example, \code{Gauge.Day.Type == "Target"}.
#' @param downwind_control_subset A logical expression used to extract the relevant subset of downwind (second stage) observations from \code{data} that were not exposed to treatment (operating ionizers). For example, \code{Gauge.Day.Type == "Control"}.
#' @param positive_subset A logical expression used to extract the relevant subset of observations from \code{data} with positive rainfall - these are the observations that are used in the fitting of upwind (first stage) LMM, downwind (second stage) LMM, downwind (second stage) treatment-only LMM, downwind (second stage) control-only LMM, and the downwind (second stage) propensity score model.
#' @param attr_type A character string specifying the type of attribution estimates. Must be one of \code{"Chambers_Chandra"}, \code{"ThoEtAl"}, or \code{"No"}. See "Details" for more information.
#' @param x_downwind_name A character vector containing variable names from the right hand side of \code{downwind_lmm_formula}, for those variables that are not related to ionizers (treatment). The intercept is always included and does not need to be specified.
#' @param target_only Logical. If \code{TRUE} the attribution estimates are computed based on only target observations. If \code{FALSE} the attribution estimates are computed based on both treatment and control observations.
#' @param bootstrap An optional logical. If \code{TRUE} bootstrap is carried out to perform inference on the attribution and sample average treatment effect. If \code{FALSE} (default) no bootstrap is carried out.
#' @param bootstrap_option An optional list containing all bootstrap settings, used only when \code{bootstrap = TRUE}. See \code{\link{bootstrap_option}} for the default list elements and their usage.
#' @param permutation An optional logical, If \code{TRUE} randomized permutation is carried out on the ionizer operation (treatment) schedule to perform inference on the attribution and sample average treatment effect. If \code{FALSE} (default) no randomized permutation is carried out.
#' @param permutation_option An optional list containing all permutation settings, used only when \code{permutation = TRUE}. See \code{\link{permutation_option}} for the default list elements and their usage.
#'
#'
#' @returns A list containing
#' \describe{
#' \item{all_fitted_models}{A list of model objects from \code{\link{lme4}{lmer}} for the first stage (upwind), second stage (downwind) LMM, second stage (downwind) treatment-only LMM, and second stage (downwind) control-only LMM, along with model objects from \code{\link{stats}{glm}} for the logistic model of rainfall event indicator (\code{NULL} if \code{downwind_logistic_formula} is not specified) and the propensity score model for the treatment indicator of second stage (downwind) observations.}
#' \item{hatattr}{A vector containing the attribution estimates.}
#' \item{hatsate}{A vector containing the sample average treatment effect estimates.}
#' \item{bootstrap_result}{A list of matrices with the following elements:}
#' \describe{
#'  \item{hatattr}{Matrix of bootstrap samples for attribution estimates.}
#'  \item{hatsate}{Matrix of bootstrap samples for SATE estimates.}
#'  \item{downwind_lmm_param}{Matrix of bootstrap samples for fixed effect coefficient and random effect variance estimates of downwind (second stage) LMM.}
#'  \item{downwind_logistic_param}{Matrix of bootstrap samples for regression coefficient estimates of downwind logistic model fitted to the rainfall event indicators. This is \code{NULL} if \code{downwind_logistic_formula} is not specified.}
#'  \item{downwind_propensity_param}{Matrix of bootstrap samples for regression coefficient estimates of downwind propensity score model fitted to the treatment indicators.}
#'  \item{downwind_positive_target_lmm_param}{Matrix of bootstrap samples for fixed effect coefficient and random effect variance estimates of downwind (second stage) treatment-only LMM.}
#'  \item{downwind_positive_control_lmm_param}{Matrix of bootstrap samples for fixed effect coefficient and random effect variance estimates of downwind (second stage) control-only LMM.}
#'  \item{downwind_LogRain}{Matrix of bootstrap samples for the log-transformed rainfall of all downwind (second-stage) observations. Observations with zero bootstrapped rainfall are represented as \code{NA}.}
#' }
#' \item{bootstrap_CI_result}{A list of matrices with same element names as in \code{bootstrap_result} (excluding \code{downwind_LogRain}), containing the corresponding bootstrap percentile confidence intervals.}
#' \item{bootstrap_p_value_result}{A list of matrices with same element names as in \code{bootstrap_result} (excluding \code{downwind_LogRain}), containing the corresponding proportion of bootstrap samples that are less than zero.}
#' \item{bootstrap_plot_result}{A list of matrices with two elements:}
#' \describe{
#' \item{hatattr}{A list of \code{ggplot} objects, each showing the bootstrap distribution of attribution estimates. Each plot includes a dotted vertical line at zero and a solid vertical line at the original estimate based on the observed data.}
#' \item{hatsate}{A list of \code{ggplot} objects, each showing the bootstrap distribution of SATE estimates. Each plot includes a dotted vertical line at zero and a solid vertical line at the original estimate based on the observed data.}
#' }
#' \item{permutation_result}{A list of matrices with the two elements:}
#' \describe{
#' \item{hatattr}{Matrix of permutation samples for attribution estimates.}
#' \item{hatsate}{Matrix of permutation samples for SATE estimates.}
#' }
#' \item{permutation_p_value_result}{A list of matrices with same element names as in \code{permutation_result}, containing the corresponding proportion of permutation samples that are greater than or equal to the original estimate based on the observed data.}
#' \item{permutation_plot_result}{A list of matrices with two elements:}
#' \describe{
#' \item{hatattr}{A list of \code{ggplot} objects, each showing the permutation distribution of attribution estimates. Each plot includes a solid vertical line at the original estimate based on the observed data.}
#' \item{hatsate}{A list of \code{ggplot} objects, each showing the permutation distribution of SATE estimates. Each plot includes a solid vertical line at the original estimate based on the observed data.}
#' }
#'}

#TODO: Fix the formatting of "Value" section, so that bootstrap_result is correctly formatted
#TODO: Change naming of attr_type to something like chambersETAL, thoETAL
rain_attr = function(data, upwind_lmm_formula, instr_pred_name, instr_pred_type,
                     downwind_lmm_formula, downwind_logistic_formula = NULL, downwind_propensity_formula,
                     rain_col_name,
                     upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset, positive_subset,
                     attr_type, x_downwind_name, target_only,
                     bootstrap = FALSE, bootstrap_option = bootstrap_option(),
                     permutation = FALSE, permutation_option = permutation_option()
                     ){

  if(!instr_pred_name %in% all.vars(downwind_lmm_formula)){
    stop("instr_pred_name cannot be found in downwind_lmm_formula")
  }

  if(formula.tools::lhs(as.formula(gsub("[()]", "", downwind_propensity_formula))) != substitute(downwind_target_subset)){
    stop("The definition of Target provided in  downwind_target_subset is not consistent with the LHS of downwind_propensity_formula")
  }

  if(!is.null(downwind_logistic_formula)){
    if(formula.tools::lhs(as.formula(gsub("[()]", "", downwind_logistic_formula))) != substitute(positive_subset)){
      stop("The definition of Positive Rainfall provided in  positive_subset is not consistent with the LHS of downwind_logistic_formula")
    }
  }else{
    if(bootstrap & bootstrap_option$bootstrap_zero){
      stop('downwind_logistic_formula needs to be specified when bootstrap = T and bootstrap_option$bootstrap_zero = T')
    }
  }


  if(mean(x_downwind_name %in% all.vars(formula.tools::rhs(downwind_lmm_formula)) )!=1 ){
    stop("At least one variable in x_downwind_name cannot be found on RHS of downwind_lmm_formula")
  }

  if(!attr_type %in% c('Chambers_Chandra','ThoEtAl', 'No')){
    stop("attr_type should be one of 'Chambers_Chandra','ThoEtAl', or 'No''")
  }

  if(bootstrap){
    if(!is.null(bootstrap_option$positive_prob_threshold) ){
      if(bootstrap_option$positive_prob_threshold < 0 | bootstrap_option$positive_prob_threshold > 1  ){
        stop("bootstrap_option$positive_prob_threshold must be either NULL or between 0 and 1")
      }
    }

    if( bootstrap_option$winsorize_individual_rain &
        (!is.numeric(bootstrap_option$individual_rain_interval) | length(bootstrap_option$individual_rain_interval)!= 2 )
        ){
      stop("bootstrap_option$individual_rain_interval must be a numeric vector of length 2")
    }

    if( bootstrap_option$winsorize_total_rain &
        (!is.numeric(bootstrap_option$total_rain_interval) | length(bootstrap_option$total_rain_interval)!= 2) ){
      stop("bootstrap_option$total_rain_interval must be a numeric vector of length 2")
    }

    if(bootstrap_option$CI_level < 0 |bootstrap_option$CI_level > 1 ){
      stop("bootstrap_option$CI_level must be between 0 and 1")
    }

  }

  if(permutation){
    if(!permutation_option$ionizer_operation_year_column_name %in% colnames(permutation_option$ionizer_operation) ){
      stop("The columns of permutation_option$ionizer_operation do not contain permutation_option$ionizer_operation_year_column_name")
    }

    if(!permutation_option$ionizer_operation_day_column_name %in% colnames(permutation_option$ionizer_operation) ){
      stop("The columns of permutation_option$ionizer_operation do not contain permutation_option$ionizer_operation_day_column_name")
    }

    if(!permutation_option$ionizer_operation_day_column_name %in% colnames(data) ){
      stop("permutation_option$ionizer_operation_day_column_name cannot be found in the column names of data")
    }

    if(mean(permutation_option$data_target_column_names %in% colnames(data)) != 1){
      stop("At least one element of permutation_option$data_target_column_names cannot be found in the column names of data")
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
  #positive = ( data[,rain_col_name] > 0)
  positive_expr = rlang::enquo(positive_subset)
  positive = rlang::eval_tidy(positive_expr, data = data)

  #Binary indicator of length N_downwind_positive, indicating whether or not each downwind positive observation is target observation
  downwind_positive_target_expr = rlang::enquo(downwind_target_subset)
  downwind_positive_target = rlang::eval_tidy(downwind_positive_target_expr, data = data[downwind & positive,])

  #Binary indicator of length N_downwind_positive, indicating whether or not each downwind positive observation is control observation
  downwind_positive_control_expr = rlang::enquo(downwind_control_subset)
  downwind_positive_control = rlang::eval_tidy(downwind_positive_control_expr, data = data[downwind & positive,])

  #Fit Upwind LMM and get the instrumental prediction, then fit Downwind LMM, and Downwind Logistic Model
  fitted_models = fit_upwind_downwind_models(data, upwind_lmm_formula, instr_pred_name, instr_pred_type, downwind_lmm_formula, downwind_logistic_formula, upwind, downwind, positive)

  downwind_positive_data = fitted_models$data[downwind & positive, ]
  #Compute Point Estimates for Attribution - using Chambers_Chandra or ThoEtAl Estimates
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
                                          individual_rain_interval = bootstrap_option$individual_rain_interval,
                                          winsorize_total_rain = bootstrap_option$winsorize_total_rain,
                                          total_rain_interval = bootstrap_option$total_rain_interval,
                                          bootstrap_seed = bootstrap_option$bootstrap_seed,
                                          bootstrap_parallel = bootstrap_option$bootstrap_parallel,
                                          bootstrap_parallel_num_worker = bootstrap_option$bootstrap_parallel_num_worker,
                                          ori_data = fitted_models$data,
                                          downwind = downwind,
                                          ori_positive = positive,
                                          rain_col_name = rain_col_name,
                                          downwind_target_expr = downwind_positive_target_expr,
                                          downwind_control_expr = downwind_positive_control_expr,
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
                                             permutation_seed = permutation_option$permutation_seed,
                                             permutation_parallel = permutation_option$permutation_parallel,
                                             permutation_parallel_num_worker = permutation_option$permutation_parallel_num_worker,
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
#Note that for attr_type == 'ThoEtAl', the hatSigma_beta matrix is ALWAYS obtained from downwind_lmm_fit
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



  if(attr_type == 'Chambers_Chandra'){
    if(is.null(hatu)){
      hatu = predict(downwind_lmm_fit, newdata = downwind_positive_data, random.only = TRUE)
    }

    #compute log_hatw = (1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t, for PDR (i,t) or PDR Target (i,t)
    #Depends on if there is offset term on LHS of formula(downwind_lmm_fit)
    if(length(formula.tools::lhs.vars(formula(downwind_lmm_fit))) == 1){
      log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row])
    }

    if(length(formula.tools::lhs.vars(formula(downwind_lmm_fit))) > 1){
      all_offset_terms = formula.tools::lhs.vars(formula(downwind_lmm_fit))[2:length(formula.tools::lhs.vars(formula(downwind_lmm_fit)))]
      if(length(all_offset_terms) > 1){
        log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row]) + apply(downwind_positive_data[,all_offset_terms],1,sum)
      }

      if(length(all_offset_terms) == 1){
        log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row]) + as.vector(downwind_positive_data[,all_offset_terms])
      }
    }

    # log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row])

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

  if(attr_type == 'ThoEtAl'){
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
  #where we replace y_i with (LogRain_i - natural_pred_i) which is captured by lme4::getME(downwind_lmm_fit, 'y') ) below that returns the response of b_downwind_lmm_fit i.e., Lograin_i^* - natural_pred_i
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



#' @title
#' Bootstrap Procedure for Rainfall Enhancement Trial Data
#'
#' @description
#' Implements the two-level bootstrap procedure used in \code{\link{rain_attr}} for inference on attribution and sample average treatment effect in rainfall enhancement trial data.
#'
#' @details
#' This function implements the bootstrap procedure used in \code{\link{rain_attr}}.
#' It is intended for internal use only. Users should not call this function
#' directly. Instead, bootstrap inference should be performed by calling \code{\link{rain_attr}} with \code{bootstrap = TRUE}.
#'
#' \strong{First-Level Bootstrap} \cr
#' The first-level is an optional level that is only carried out when \code{bootstrap_zero = TRUE}. This level generates bootstrap samples of binary rainfall event indicator \eqn{L_{ij}^*} via \eqn{P(L_{ij}^* = 1) = \hat{\gamma}_{ij} } for the subset of observations from \code{ori_data} satisfying \code{downwind}, where \eqn{\hat{\gamma}_{ij}} are the predicted probabilities from the downwind logistic model fitted to the original binary rainfall event indicators \eqn{L_{ij}} i.e., \eqn{\hat{\gamma}_{ij}} = \code{predict(ori_fitted_models$downwind_logistic_fit,type = "response")}.
#' It is worth noting that if \code{positive_prob_threshold} is supplied, then \eqn{\hat{\gamma}_{ij}} that are less than \code{positive_prob_threshold} are set to be zeros before being used to generate \eqn{L_{ij}^*}.
#' When \code{bootstrap_zero = FALSE}, the first-level bootstrap is not carried out and thus \eqn{L_{ij}^* = L_{ij}}.
#'
#' \strong{Second-Level Bootstrap} \cr
#' Recall that the original downwind (second stage) LMM is given as
#' \deqn{y_{ij} = x_{ij}^\top \alpha + z_{ij}^\top \beta + u_i + e_{ij},}
#' the second-level bootstrap generates bootstrap samples of positive rainfall for the subset of observations from \code{ori_data} satisfying \code{downwind} and \eqn{L_{ij}^* = 1} from the above first-level bootstrap via
#' \deqn{y_{ij}^* = x_{ij}^\top \hat{\alpha} + z_{ij}^\top \hat{\beta} + u_i^* + e_{ij}^* ,}
#' where \eqn{\hat{\alpha}} and \eqn{\hat{\beta}} are estimated fixed effect coefficients from the fitted downwind (second stage) LMM i.e., \code{ori_fitted_models$downwind_lmm_fit}, and \eqn{u_i^*} and \eqn{e_{ij}^*} are bootstrap samples of the random intercepts and error terms.
#'
#' We now explain the procedure for obtaining bootstrap samples of \eqn{u_i^*} and \eqn{e_{ij}^*}.
#' Let \eqn{\hat{r}_{ij} = y_{ij} - x_{ij}^\top \hat{\alpha} - z_{ij}^\top \hat{\beta}} be the marginal residuals for the subset of observations from \code{ori_data} satisfying \code{downwind & ori_positive}, \eqn{\hat{u}_i = \sum_{j=1}^{n_i} \hat{r}_{ij}/n_i } be the day(group)-level average residuals, and \eqn{\hat{e}_{ij} = r_{ij} - \hat{u}_i } be the gauge(unit)-level residuals, where \eqn{n_i} denotes the total number of observations in day (group) \eqn{i} from \code{ori_data} satisfying \code{downwind & ori_positive}.
#' Depending on the chosen \code{bootstrap_type}, bootstrap samples of \eqn{u_i^*} and \eqn{e_{ij}^*} are generated via:
#' \describe{
#' \item{\code{REB0}}{
#'    \itemize{
#'     \item \eqn{u_i^* = SRSWR( (\hat{u}_1, \ldots, \hat{u}_D ), 1 )} for \eqn{i=1,\ldots, D^*}, where \eqn{D} is the number of unique days (groups) from \code{ori_data} satisfying \code{downwind & ori_positive}, \eqn{D^*} is the number of unique days (groups) from \code{ori_data} satisfying \code{downwind} and \eqn{L_{ij}^* = 1}, and \eqn{SRSWR(a,c)} denote the outcome of \eqn{c} independent draws based on simple random sampling with replacement from the vector \eqn{a}.
#'     \item First, sample the donor cluster \eqn{d_i^* = SRSWR( ( 1,\ldots,D ), 1 ) } for \eqn{i = 1,\ldots, D^*}. Then, sample
#'     \eqn{e_i^* = (e_{i1}^*, \ldots, e_{in_i^*}^*)^\top = SRSWR( ( \hat{e}_{d_i^* 1}, \ldots, \hat{e}_{d_i^* n_{d_i^*}} ), n_i^*    )  }, where \eqn{n_i^*} denotes the total number of observations in day (group) \eqn{i} from \code{ori_data} satisfying \code{downwind} and \eqn{L_{ij}^* = 1}.
#'    }
#' }
#'
#' \item{\code{REB1}}{
#' Replaces \eqn{\hat{u}_i} and \eqn{\hat{e}_{ij}} in \code{REB0} with
#' \eqn{\hat{u}_{ij}^{cs} = \hat{\sigma}_u \hat{u}_i^{c} \{ D^{-1} \sum_{i' =1}^{D} \hat{u}_i^2 \}^{-1/2} } and
#' \eqn{\hat{e}_{ij}^{s} = \hat{\sigma}_e \hat{e}_{ij} \{ N^{-1} \sum_{i' = 1}^{D} \sum_{j' = 1}^{n_{i'}  } \hat{e}_{i'j'}^2  \}^{-1/2} }, respectively,
#' where \eqn{ \hat{u}_i^c = \hat{u}_i - D^{-1} \sum_{i'=1}^{D} \hat{u}_{i'} }, \eqn{\hat{\sigma}^2_u} and \eqn{\hat{\sigma}^2_e} are the estimated variances of random intercepts and error terms from the fitted downwind (second stage) LMM, and \eqn{N = \sum_{i=1}^{D} n_i} is the total number of observations from \code{ori_data} satisfying \code{downwind & ori_positive}.
#' }
#'
#' \item{\code{REB2}}{
#' Same procedure for obtaining \eqn{u_i^*} and \eqn{e_{ij}^*} as in \code{REB0}, but involves an additional post-processing step on the bootstrap estimates of attribution and SATE, as discussed below.
#' }
#'
#' \item{\code{PREB0}}{
#'    \itemize{
#'     \item \eqn{u_i^* = SRSWR( (\hat{u}_1, \ldots, \hat{u}_D ), 1 )} for \eqn{i=1,\ldots, D^*}.
#'     \item First, sample the donor cluster \eqn{d_i^* = PPSWR{ (1,\ldots,D), (n_1,\cdots, n_D), 1 } } for \eqn{i = 1,\ldots, D^*},
#'     where \eqn{PPSWR(a,b,c)} denotes the outcome of \eqn{c} independent draws based on probability-proportional-to-size sampling with replacement from the vector \eqn{a = (a_1,\ldots, a_D)} with corresponding sizes given by the vector \eqn{b = (b_1,\ldots,b_D)},
#'     i.e., the probability of \eqn{a_i} being selected is given as \eqn{b_i / \sum_{i' = 1}^{D} b_{i'}}.
#'     Then, sample
#'     \eqn{e_i^* = (e_{i1}^*, \ldots, e_{in_i^*}^*)^\top = SRSWR( ( \hat{e}_{d_i^* 1}, \ldots, \hat{e}_{d_i^* n_{d_i^*}} ), n_i^*    )  }.
#'    }
#' }
#'
#' \item{\code{PREB1}}{
#' Replaces \eqn{\hat{u}_i} and \eqn{\hat{e}_{ij}} in \code{PREB0} with
#' \eqn{\hat{u}_{ij}^{sc} = \hat{\sigma}_u \hat{u}_i^{c} \{ D^{-1} \sum_{i' =1}^{D} (\hat{u}_i^c)^2 \}^{-1/2} } and
#' \eqn{\hat{e}_{ij}^{s} = \hat{\sigma}_e \hat{e}_{ij} \{ N^{-1} \sum_{i' = 1}^{D} \sum_{j' = 1}^{n_{i'}  } \hat{e}_{i'j'}^2  \}^{-1/2} }, respectively.
#' }
#'
#' \item{\code{PREB2}}{
#' Same procedure for obtaining \eqn{u_i^*} and \eqn{e_{ij}^*} as in \code{PREB0}, but involves an additional post-processing step on the bootstrap estimates of attribution and SATE, as discussed below.
#' }
#'
#' \item{\code{MREB1}}{
#' Replaces \eqn{\hat{u}_i} and \eqn{\hat{e}_{ij}} in \code{REB0} with
#' \eqn{\hat{u}_{ij}^{sc}} defined under \code{PREB1} and
#' \eqn{\tilde{e}_{ij}^{s} = \hat{\sigma}_e \hat{e}_{ij} \{  \sum_{i' = 1}^{D} \sum_{j' = 1}^{n_{i'}  } D^{-1} n_{i'}^{-1} \hat{e}_{i'j'}^2  \}^{-1/2} }, respectively.
#' }
#'
#' }
#'
#' The random effect block (REB0, REB1, REB2) bootstraps proposed in Chambers and Chandra (2013) were originally designed to handle balanced clustered data, while the proportional REB (PREB0, PREB1, PREB2) bootstraps and the MREB1 bootstrap proposed by Tho et al. (2025) are generalizations of the REB bootstraps to accommodate highly unbalanced clustered data.
#' Therefore, it is recommended to use either \code{bootstrap_type = "PREB1"} or \code{bootstrap_type = "MREB1"}, especially when \eqn{n_i}'s are highly unbalanced. Users are refered to Tho et al. (2025) for more discussion on the comparison among these bootstrap methods.
#'
#' \strong{Adjustment of Bootstrapped Rainfall} \cr
#' After obtaining the bootstrapped \eqn{y_{ij}^*} that are assumed to be log-rainfall, it is possible to perform some adjustment to the bootstrap samples to ensure that they are similar to the observed data.
#'
#' Let \eqn{Rain_{ij}^* = \exp(y_{ij}^*)} be the bootstrapped raw rainfall.
#' When \code{discretize_rain = TRUE}, bootstrap samples that satisfy \eqn{ Rain_{ij}^* \in (0, 0.3] } are replaced as 0.2, \eqn{ Rain_{ij}^* \in (0.3, 0.5] } are replaced as 0.4, \eqn{ Rain_{ij}^* \in (0.5, 0.7] } are replaced as 0.6, and \eqn{ Rain_{ij}^* \in (0.7, 0.9] } are replaced as 0.8.
#'
#' When \code{winsorize_individual_rain = TRUE}, bootstrap samples that satisfy \eqn{Rain_{ij}^* >  } \code{individual_rain_upper} are replaced as \code{individual_rain_center} + ( \code{individual_rain_upper} -  \code{individual_rain_center} ) \eqn{\times U_{ij}^* }, where each of the \eqn{U_{ij}^*} are i.i.d. standard uniform random numbers.
#'
#' When \code{winsorize_total_rain = TRUE}, if \eqn{ \sum_{(i,j)} Rain_{ij}^* \notin [} \code{total_rain_lower}, \code{total_rain_upper} \eqn{]} where the summation is over the subset of observations in \code{ori_data} satisfying \code{downwind} and \eqn{L_{ij}^* = 1}, then all bootstrap samples of \eqn{Rain_{ij}^*} are rescaled by a common factor of \code{runif(n = 1, min = total_rain_lower, max = total_rain_upper)} \eqn{ / \sum_{(i,j)} Rain_{ij}^*  }.
#'
#' The final adjusted \eqn{Rain_{ij}^*} are converted back to the log-scale based on the formula \eqn{y_{ij}^* = \log(Rain_{ij}^*)}. Therefore, these adjustments should only be used when modelling log-transformed rainfall in the two-stage LMM approach, but not raw rainfall.
#' If an offset term is included on the LHS of \code{downwind_lmm_formula} e.g., \code{downwind_lmm_formula = LogRain - Offset}, the above adjustments could still be used as the function would use the relationship \eqn{Rain_{ij}^* = \exp( y_{ij}^* + Offset_{ij} )} and \eqn{y_{ij}^* = \log(Rain_{ij}^*) + Offset_{ij} }.
#'
#' \strong{Bootstrap Distributions of Attribution, SATE and other parameters} \cr
#' The same procedure described in \code{\link{rain_attr}} is then repeated on the bootstrapped dataset, where the bootstrap samples of \eqn{y_{ij}^*} generated from the above two-level bootstrap procedure are used to replace the original observed \eqn{y_{ij}}.
#' This includes the fitting of the downwind (second stage) LMM, downwind (second-stage) target-only LMM, and downwind (second-stage) control-only LMM to the subset of observations from \code{ori_data} satisfying \code{downwind &} \eqn{ L_{ij}^* = 1}. When \code{bootstrap_zero = TRUE}, two additional models are also fitted; namely,
#' the downwind logistic model fitted to the subset of observations from \code{ori_data} satisfying \code{downwind} with the response being the bootstrapped rainfall event indicator \eqn{L_{ij}^*}, and
#' the downwind propensity score model fitted to the subset of observations from \code{ori_data} satisfying \code{downwind &} \eqn{L_{ij}^* = 1} with the response being the indicator \eqn{I_{ij}} for exposure to ionizer (treatment).
#' Finally, two attribution estimates and five SATE estimates are computed based on the estimation results of these models fitted to the bootstrapped dataset, where \eqn{Rain_{ij}} and \eqn{y_{ij}} are replaced by their bootstrap counterparts \eqn{Rain_{ij}^*} and \eqn{y_{ij}^*}, respectively.
#'
#' By repeatedly generating bootstrap samples, fitting models and computing attribution and SATE estimates for \code{B_bootstrap} number of times, this function returns the bootstrap distributions of
#' \itemize{
#'    \item{Two attribution estimates}
#'    \item{Five SATE estimates}
#'    \item{Fixed effect coefficients and random effect variance estimates of downwind LMM, downwind treatment-only LMM, and downwind control-only LMM.}
#'    \item{Log-transformed rainfall \eqn{\log(Rain_{ij}^*)}}
#'    \item{Regression coefficient estimates of downwind logistic model fitted to the bootstrapped rainfall event indicators, only when \code{bootstrap_zero = TRUE}}
#'    \item{Regression coefficient estimates of downwind propensity score model, only when \code{bootstrap_zero = TRUE}}
#' }
#'
#' \strong{Parallel Bootstrap Execution} \cr
#' This function also supports parallel execution via the argument \code{bootstrap_parallel}.  When \code{bootstrap_parallel = TRUE}, the \code{B_bootstrap} repeated iterations are executed
#' concurrently across \code{bootstrap_parallel_num_worker} workers. Each iteration, including rainfall resampling, model fitting, and parameter estimation, is performed independently on different workers, which can reduce computation time compared to sequential execution.
#' The parallel execution is implemented using the \code{\link[foreach]{foreach}} package with the \code{\link[doParallel]{doParallel}} backend and \code{\link[doRNG]{registerDoRNG}} for reproducibility.
#'
#' \strong{Post-processing of \code{REB2} and \code{PREB2}} \cr
#' The \code{REB2} and \code{PREB2} perform post-processing on all of the above bootstrap distributions
#' (except for log-transformed rainfall) generated by \code{REB0} and \code{PREB0}, respectively.
#'
#' Specifically:
#' \itemize{
#'   \item The bootstrap distributions of attribution, SATE, fixed-effect coefficients from LMMs,
#'   and regression coefficients from the downwind logistic model and downwind propensity score model
#'   are mean-centered at their original estimates using:
#'   \deqn{\hat{\theta}^*_b \leftarrow \hat{\theta}^*_b + \hat{\theta} - \frac{1}{B}\sum_{b=1}^{B}\hat{\theta}^*_b}
#'   where \eqn{\hat{\theta}^*_b} is the estimate from the \eqn{b}-th bootstrap sample and
#'   \eqn{\hat{\theta}} is the corresponding estimate from the original dataset.
#'
#'   \item For LMM random effect variance components, the bootstrap distributions of random intercept
#'   variance (\eqn{\sigma^2_u}) and residual variance (\eqn{\sigma^2_e}) are first adjusted to be
#'   empirically uncorrelated (see Section 2.3.3 of Chambers & Chandra (2013), and Section 2.3 of Tho et al. (2025)).
#'
#'   \item After adjustment, ratio corrections are applied:
#'   \deqn{\hat{\sigma}^{2*}_{u,b} \leftarrow \hat{\sigma}^{2*}_{u,b} \times \frac{\hat{\sigma}^2_u}{\frac{1}{B}\sum_{b=1}^{B}\hat{\sigma}^{2*}_{u,b}}}
#'   \deqn{\hat{\sigma}^{2*}_{e,b} \leftarrow \hat{\sigma}^{2*}_{e,b} \times \frac{\hat{\sigma}^2_e}{\frac{1}{B}\sum_{b=1}^{B}\hat{\sigma}^{2*}_{e,b}}}
#'   where \eqn{\hat{\sigma}^{2*}_{u,b}} and \eqn{\hat{\sigma}^{2*}_{e,b}} are the adjusted bootstrap estimates,
#'   and \eqn{\hat{\sigma}^2_u}, \eqn{\hat{\sigma}^2_e} are the original estimates.
#' }
#'
#' \strong{Bootstrap P-values, Bootstrap Percentile Confidence Intervals, Bootstrap Plots} \cr
#' All bootstrap distributions (except for log-transformed rainfall) produced by this function are further used in \code{\link{rain_attr}} to:
#' \itemize{
#'
#' \item{Compute bootstrap p-values as the proportion of bootstrapped estimates that are less than zero,
#' i.e., \deqn{ \frac{1}{B} \sum_{b=1}^{B} 1_{ \{ \hat{\theta}_b^* < 0 \}} }
#' where \eqn{1_{\{\cdot\}}} is the indicator function and \eqn{\hat{\theta}_b} denotes the estimate from the \eqn{b}-th bootstrap sample.
#'  }
#'
#'  \item{Compute bootstrap percentile confidence interval with confidence level \eqn{(1-\alpha) \times 100\%} as
#'  \deqn{ [ \hat{\theta}^*_{\alpha / 2}, \hat{\theta}^*_{1 - \alpha/2} ], }
#'  where \eqn{\hat{\theta}^*_{p}} is the empirical \eqn{p}-th quantile of the bootstrap distribution \eqn{\{ \hat{\theta}^*_{1}, \ldots, \hat{\theta}^*_{B} \}} of the parameter estimate \eqn{\hat{\theta}}
#'   }
#' }
#'
#' Additionally, the bootstrap distributions of attribution and SATE are also plotted using kernel density estimates, with:
#' \itemize{
#'   \item A solid vertical line for the original estimate \eqn{\hat{\theta}};
#'   \item A dotted vertical line at zero for reference.
#' }
#'
#' @param B_bootstrap An integer specifying the number of bootstrap replicates.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param bootstrap_type A character string specifying the type of bootstrap.
#'   Must be one of \code{"REB0"}, \code{"REB1"}, \code{"REB2"}, \code{"PREB0"},
#'   \code{"PREB1"}, \code{"PREB2"}, or \code{"MREB1"}.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param bootstrap_zero Logical. If \code{TRUE}, the optional first-level bootstrap is performed to generate bootstrap samples of binary rainfall event indicators.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param positive_prob_threshold An optional numeric value between 0 and 1 specifying the probability threshold for generating bootstrap samples of binary rainfall event indicators. Probabilities below this threshold are set to zero.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param discretize_rain Logical. If \code{TRUE}, rainfall values are discretized in bootstrap resamples.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param winsorize_individual_rain Logical. If \code{TRUE}, individual rainfall values in bootstrap samples that exceed the upper bound specified by \code{individual_rain_interval} are replaced with random draws from a uniform distribution over \code{[individual_rain_interval[1], individual_rain_interval[2]]}.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param individual_rain_interval Numeric vector of length 2 specifying the lower and upper bounds for adjusting bootstrapped individual rainfall values that are too large when \code{winsorize_individual_rain = TRUE}.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param winsorize_total_rain Logical. If \code{TRUE}, all individual rainfall values in each bootstrap sample are proportionally rescaled so that the total equals a random number drawn uniformly from
#'   \code{[total_rain_interval[1], total_rain_interval[2]]} whenever the total bootstrapped rainfall falls outside this interval.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param total_rain_interval Numeric vector of length 2 specifying the lower and upper bounds for adjusting the total of bootstrapped rainfall values when \code{winsorize_total_rain = TRUE}.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param bootstrap_seed n integer specifying the random seed for the bootstrap procedure. Reproducibility is guaranteed only if \code{bootstrap_parallel} is the same, since parallel execution changes the order of random number generation.
#' (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param bootstrap_parallel Logical. If \code{TRUE}, each bootstrap run is executed in parallel across multiple workers. If \code{FALSE}, they are run sequentially.
#' (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#' @param bootstrap_parallel_num_worker An integer specifying the number of parallel workers to use when \code{bootstrap_parallel = TRUE}.
#' (User-configurable bootstrap option using \code{\link{bootstrap_option}})
#'
#' @param ori_data A data frame containing the original dataset used in \code{\link{rain_attr}}, along with an additional column containing the fitted values generated from the upwind (first stage) LMM.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind A logical vector indicating which observation in \code{ori_data} would be used in the downwind (second stage) LMM fitting.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param ori_positive A logical vector indicating which observation in \code{ori_data} has positive rainfall.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param rain_col_name A character string specifying the column name of the raw scale rainfall in \code{ori_data}.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind_target_expr A quosure (created using `rlang::enquo()`) representing a logical expression used to extract the relevant subset of downwind (second stage) observations from \code{ori_data} that were exposed to treatment (operating ionizers).
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind_control_expr A quosure (created using `rlang::enquo()`) representing a logical expression used to extract the relevant subset of downwind (second stage) observations from \code{ori_data} that were not exposed to treatment (operating ionizers).
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param ori_fitted_models A list containing the models fitted to the \code{ori_data}, including the upwind (first stage) LMM, downwind (second stage) LMM, downwind (second stage) treatment-only LMM, downwind (second stage) target-only LMM, downwind (second stage) logistic model for rainfall event indicator, and downwind (second stage) propensity score model.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind_lmm_formula A two sided linear formula object to be used in \link[lme4]{lmer}, describing both the fixed-effects and random intercept part of the downwind (second stage) LMM.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param attr_type A character string specifying the type of attribution estimates. Must be one of \code{"Chambers_Chandra"}, \code{"ThoEtAl"}, or \code{"No"}. See \code{\link{rain_attr}} for more information.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param x_downwind_name A character vector containing variable names from the right hand side of \code{downwind_lmm_formula}, for those variables that are not related to ionizers (treatment). The intercept is always included and does not need to be specified.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param target_only Logical. If \code{TRUE} the attribution estimates are computed based on only target observations. If \code{FALSE} the attribution estimates are computed based on both treatment and control observations.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind_propensity_formula A two sided linear formula object to be used in \code{\link{stats}{glm}} with \code{family = "binomial"}, for fitting a propensity score model to the treatment indicators of downwind (second stage) observations.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param ori_attr_est A numeric vector containing the original attribution estimates (\code{apo} and \code{apl}) from the original dataset.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param ori_sate_est A numeric vector containing the original SATE estimates (\code{sate.mb}, \code{sate.ipw}, \code{sate.ipw.l}, \code{sate.ipw.ma} and \code{sate.aipw}) from the original dataset.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})

#' @returns A list containing
#' \describe{
#'  \item{hatattr}{Matrix of bootstrap samples for attribution estimates.}
#'  \item{hatsate}{Matrix of bootstrap samples for SATE estimates.}
#'  \item{downwind_lmm_param}{Matrix of bootstrap samples for fixed effect coefficient and random effect variance estimates of downwind (second stage) LMM.}
#'  \item{downwind_logistic_param}{Matrix of bootstrap samples for regression coefficient estimates of downwind logistic model fitted to the rainfall event indicators. This is \code{NULL} if \code{bootstrap_zero = FALSE}.}
#'  \item{downwind_propensity_param}{Matrix of bootstrap samples for regression coefficient estimates of downwind propensity score model fitted to the treatment indicators. This is \code{NULL} if \code{bootstrap_zero = FALSE}.}
#'  \item{downwind_positive_target_lmm_param}{Matrix of bootstrap samples for fixed effect coefficient and random effect variance estimates of downwind (second stage) treatment-only LMM.}
#'  \item{downwind_positive_control_lmm_param}{Matrix of bootstrap samples for fixed effect coefficient and random effect variance estimates of downwind (second stage) control-only LMM.}
#'  \item{downwind_LogRain}{Matrix of bootstrap samples for the log-transformed rainfall of all downwind (second-stage) observations. Observations with zero bootstrapped rainfall are represented as \code{NA}.}
#' }
#'
#' @seealso \code{\link{rain_attr}} for the main function, \code{\link{bootstrap_option}} for specifying bootstrap options

#TODO: Consider to add another variation of 'PREB2' and 'REB2' for adjusting downwind_lmm_fit's fixef as well as random effects, and plug these corrected estimates to compute hatattr and hatsate, rather than directly centering hatattr and hatsate
#TODO: Consider to not having to refit downwind_propernsity_formula in every bootstrap run if bootstrap_zero = FALSE, since they should be exactly the same when bootstrap_zero = FALSE
#TODO: Think if we really need to save the bootstrapped LogRain for each bootstrap run, as this would be a LARGE object when B_bootstrap is large since it is of dimension B_bootstrap x num_downwind (43276)
bootstrap_downwind = function(B_bootstrap, bootstrap_type, bootstrap_zero, positive_prob_threshold = NULL, discretize_rain, winsorize_individual_rain, individual_rain_interval, winsorize_total_rain, total_rain_interval,
                              bootstrap_seed, bootstrap_parallel, bootstrap_parallel_num_worker,
                              ori_data, downwind, ori_positive, rain_col_name, downwind_target_expr, downwind_control_expr, ori_fitted_models,
                              downwind_lmm_formula, attr_type, x_downwind_name, target_only,
                              downwind_propensity_formula,
                              ori_attr_est, ori_sate_est){

  if(!bootstrap_type %in% c('REB0','REB1','REB2','PREB0','PREB1','PREB2','MREB1')){
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
  if(bootstrap_zero){
    bootstrap_downwind_logistic_param_matrix = matrix(data = NA, nrow = B_bootstrap,
                                                                       ncol = length(coef(ori_fitted_models$downwind_logistic_fit)),
                                                                       dimnames = list(NULL, names(coef(ori_fitted_models$downwind_logistic_fit))))

    bootstrap_downwind_propensity_param_matrix = matrix(data = NA, nrow = B_bootstrap,
                                                        ncol = length(coef(ori_fitted_models$downwind_propensity_fit)),
                                                        dimnames = list(NULL, names(coef(ori_fitted_models$downwind_propensity_fit))))
    }

  bootstrap_downwind_positive_target_lmm_param_matrix = matrix(data = NA, nrow = B_bootstrap,
                                                               ncol = length(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'vcov']),
                                                               dimnames = list(NULL, c(names(lme4::fixef(ori_fitted_models$downwind_positive_target_lmm_fit)), paste0('VarComponent_',as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_target_lmm_fit))[,'grp'] ) )))
  bootstrap_downwind_positive_control_lmm_param_matrix = matrix(data = NA, nrow = B_bootstrap,
                                                                ncol = length(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)) + length(as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'vcov']),
                                                                dimnames = list(NULL, c(names(lme4::fixef(ori_fitted_models$downwind_positive_control_lmm_fit)), paste0('VarComponent_',as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_positive_control_lmm_fit))[,'grp'] ) )))

  bootstrap_downwind_LogRain_matrix = matrix(data = NA, nrow = B_bootstrap, ncol = num_downwind)

  z_downwind_name = setdiff(names(lme4::fixef(ori_fitted_models$downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
  downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)

  #browser()
  if(!bootstrap_parallel){
    set.seed(bootstrap_seed)
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

        #browser()

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

        #Back-transform into raw rainfall, depending on whether offset term is included on LHS of downwind_lmm_formula
        if(length(formula.tools::lhs.vars(downwind_lmm_formula)) == 1){
          b_raw_y = exp(b_y)
        }

        if(length(formula.tools::lhs.vars(downwind_lmm_formula)) > 1){
          all_offset_terms = formula.tools::lhs.vars(downwind_lmm_formula)[2:length(formula.tools::lhs.vars(downwind_lmm_formula))]
          if(length(all_offset_terms) > 1){
            b_raw_y = exp(b_y + apply(b_downwind_positive_data[,all_offset_terms],1,sum))
          }

          if(length(all_offset_terms) == 1){
            b_raw_y = exp(b_y + as.vector(b_downwind_positive_data[,all_offset_terms]))
          }
        }

        # b_raw_y = exp(b_y)


        #Perform (optional) adjustment of raw rainfall
        if(discretize_rain){
          b_raw_y[b_raw_y<0.3] <- 0.2
          b_raw_y[(b_raw_y>0.3)&(b_raw_y<0.5)] <- 0.4
          b_raw_y[(b_raw_y>0.5)&(b_raw_y<0.7)] <- 0.6
          b_raw_y[(b_raw_y>0.7)&(b_raw_y<0.9)] <- 0.8
        }



        if(winsorize_individual_rain){
          b_raw_y[b_raw_y> individual_rain_interval[2]] <- individual_rain_interval[1] + (individual_rain_interval[2] - individual_rain_interval[1]) * runif(n=sum(b_raw_y> individual_rain_interval[2]))
        }

        if(winsorize_total_rain){
          if(sum(b_raw_y)< total_rain_interval[1] | sum(b_raw_y)> total_rain_interval[2]){
            b_raw_y <- b_raw_y*(runif(n=1,min=total_rain_interval[1], max=total_rain_interval[2]))/sum(b_raw_y)
          }
        }

        b_downwind_positive_data[,rain_col_name] = b_raw_y


        #Replace the original LogRain column with the newly bootstrapped (and potentially adjusted) LogRain^*
        b_downwind_positive_data[,formula.tools::lhs.vars(downwind_lmm_formula)[1]] = log(b_raw_y)



        # b_y = log(b_raw_y)

        # #This part reconstructs the bootstrapped raw rain, depending on whether there are any offset terms specified in downwind_lmm_formula
        # if(length(formula.tools::lhs.vars(downwind_lmm_formula)) == 1){
        #   b_downwind_positive_data[,rain_col_name] = exp(b_y)
        # }
        #
        # if(length(formula.tools::lhs.vars(downwind_lmm_formula)) > 1){
        #   all_offset_terms = formula.tools::lhs.vars(downwind_lmm_formula)[2:length(formula.tools::lhs.vars(downwind_lmm_formula))]
        #   if(length(all_offset_terms) > 1){
        #     b_downwind_positive_data[,rain_col_name] = exp(b_y + apply(b_downwind_positive_data[,all_offset_terms],1,sum))
        #   }
        #
        #   if(length(all_offset_terms) == 1){
        #     b_downwind_positive_data[,rain_col_name] = exp(b_y + as.vector(b_downwind_positive_data[,all_offset_terms]))
        #   }
        # }



        #Create a new column 'bootstrapped_y' instead of replacing the LogRain column to accommodate for the case of having LogRain - natural_pred on the LHS of downwind_lmm_formula
        #In this case, we are essentially creating bootstrapped_y = LogRain* - natural_pred, where LogRain* = natural_pred + Xhatbeta + u* + e*
        # b_downwind_positive_data$bootstrapped_y = b_y
        # b_downwind_lmm_fit = lme4::lmer(update.formula(downwind_lmm_formula, bootstrapped_y ~ . ),
        #                                 data = b_downwind_positive_data)

        b_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula,
                                        data = b_downwind_positive_data)

        bootstrap_downwind_lmm_param_matrix[b,] = c(lme4::fixef(b_downwind_lmm_fit),
                                                    as.data.frame(lme4::VarCorr(b_downwind_lmm_fit))[,'vcov'])



        bootstrap_downwind_LogRain_matrix[b, b_downwind_positive] = log(b_raw_y)
        bootstrap_downwind_LogRain_matrix[b, !b_downwind_positive] = NA


        # b_downwind_positive_target = b_downwind_positive_data$Gauge.Day.Type == 'Target'
        # b_downwind_positive_control = b_downwind_positive_data$Gauge.Day.Type == 'Control'
        #browser()
        b_downwind_positive_target = rlang::eval_tidy(downwind_target_expr, data = b_downwind_positive_data)
        b_downwind_positive_control = rlang::eval_tidy(downwind_control_expr, data = b_downwind_positive_data)


        b_hatattr = attr_est(attr_type, b_downwind_positive_data, rain_col_name, b_downwind_positive_target, b_downwind_positive_control,
                             x_downwind_name, target_only, downwind_lmm_fit = b_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

        # #Need to also update the downwind_separate_formula here, since the bootstrapped 2nd stage response is now stored in the column 'bootstrapped_y'
        # b_hatsate = sate_est(b_downwind_positive_data, b_downwind_positive_target, b_downwind_positive_control, downwind_propensity_formula, update.formula(downwind_separate_formula, bootstrapped_y ~ . ),
        #                      x_downwind_name, downwind_lmm_fit = b_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

        b_hatsate = sate_est(b_downwind_positive_data, b_downwind_positive_target, b_downwind_positive_control, downwind_propensity_formula, downwind_separate_formula,
                             x_downwind_name, downwind_lmm_fit = b_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

        if(bootstrap_zero){bootstrap_downwind_propensity_param_matrix[b,] = coef(b_hatsate$fitted_models$downwind_propensity_fit)}
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
  }else{
    cl = parallel::makeCluster(bootstrap_parallel_num_worker)
    #No need to use parallel:clusterExport as it does not result in any significant computational gain, nor do using foreach::foreach(..., .exports = ...)
    # if(bootstrap_export_all){
    #   parallel::clusterExport(cl, varlist = c('bootstrap_zero', 'num_downwind', 'downwind_positive_prob', 'ori_fitted_models', 'ori_data', 'downwind', 'ori_positive',
    #                                           'ori_downwind_data', 'group_name', 'ori_downwind_positive_group_label', 'cluster_sample_prob',
    #                                           'final.hat.u', 'final.hat.e', 'downwind_lmm_formula', 'rain_col_name','discretize_rain', 'winsorize_individual_rain',
    #                                           'winsorize_total_rain', 'downwind_target_expr', 'downwind_control_expr',
    #                                           'attr_type', 'x_downwind_name', 'target_only', 'downwind_propensity_formula', 'downwind_separate_formula'),
    #                           envir = environment())
    # }
    parallel::clusterExport(cl, varlist = c('attr_est', 'sate_est'))
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(seed = bootstrap_seed)
    `%dopar%` <- foreach::`%dopar%`

    results = foreach::foreach(b = 1:B_bootstrap, .packages = c("lme4","rlang","formula.tools"), .errorhandling = 'remove') %dopar% {
      if(bootstrap_zero){
        b_downwind_positive = (runif(num_downwind, min = 0, max = 1) < downwind_positive_prob)
        b_downwind_logistic_fit = glm(b_downwind_positive ~ model.matrix(ori_fitted_models$downwind_logistic_fit$formula, data = ori_data[downwind,]) - 1, data = ori_data[downwind,], family = 'binomial')
        downwind_logistic_param_b = coef(b_downwind_logistic_fit)
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

      #Back-transform into raw rainfall, depending on whether offset term is included on LHS of downwind_lmm_formula
      if(length(formula.tools::lhs.vars(downwind_lmm_formula)) == 1){
        b_raw_y = exp(b_y)
      }

      if(length(formula.tools::lhs.vars(downwind_lmm_formula)) > 1){
        all_offset_terms = formula.tools::lhs.vars(downwind_lmm_formula)[2:length(formula.tools::lhs.vars(downwind_lmm_formula))]
        if(length(all_offset_terms) > 1){
          b_raw_y = exp(b_y + apply(b_downwind_positive_data[,all_offset_terms],1,sum))
        }

        if(length(all_offset_terms) == 1){
          b_raw_y = exp(b_y + as.vector(b_downwind_positive_data[,all_offset_terms]))
        }
      }



      #Perform (optional) adjustment of raw rainfall
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

      b_downwind_positive_data[,rain_col_name] = b_raw_y


      #Replace the original LogRain column with the newly bootstrapped (and potentially adjusted) LogRain^*
      b_downwind_positive_data[,formula.tools::lhs.vars(downwind_lmm_formula)[1]] = log(b_raw_y)

      b_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula,
                                      data = b_downwind_positive_data)

      downwind_lmm_param_b = c(lme4::fixef(b_downwind_lmm_fit),
                               as.data.frame(lme4::VarCorr(b_downwind_lmm_fit))[,'vcov'])


      downwind_LogRain_b = rep(NA, num_downwind)
      downwind_LogRain_b[b_downwind_positive] = log(b_raw_y)


      b_downwind_positive_target = rlang::eval_tidy(downwind_target_expr, data = b_downwind_positive_data)
      b_downwind_positive_control = rlang::eval_tidy(downwind_control_expr, data = b_downwind_positive_data)


      b_hatattr = attr_est(attr_type, b_downwind_positive_data, rain_col_name, b_downwind_positive_target, b_downwind_positive_control,
                           x_downwind_name, target_only, downwind_lmm_fit = b_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)


      b_hatsate = sate_est(b_downwind_positive_data, b_downwind_positive_target, b_downwind_positive_control, downwind_propensity_formula, downwind_separate_formula,
                           x_downwind_name, downwind_lmm_fit = b_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

      if(bootstrap_zero){downwind_propensity_param_b = coef(b_hatsate$fitted_models$downwind_propensity_fit)}
      downwind_positive_target_lmm_param_b = c(lme4::fixef(b_hatsate$fitted_models$downwind_positive_target_lmm_fit),
                                               as.data.frame(lme4::VarCorr(b_hatsate$fitted_models$downwind_positive_target_lmm_fit))[,'vcov'])
      downwind_positive_control_lmm_param_b = c(lme4::fixef(b_hatsate$fitted_models$downwind_positive_control_lmm_fit),
                                                as.data.frame(lme4::VarCorr(b_hatsate$fitted_models$downwind_positive_control_lmm_fit))[,'vcov'])


      #TODO: try to save tor and tlr as well, to understand why the bootstrap plots for apo and apl are always the same shape
      #Maybe can also go back and look at previous plots to see if we always have same shape for apo and apl
      attr_b = c(b_hatattr$apo, b_hatattr$apl)
      sate_b = c(b_hatsate$estimates$sate.mb, b_hatsate$estimates$sate.ipw, b_hatsate$estimates$sate.ipw.l, b_hatsate$estimates$sate.ipw.ma, b_hatsate$estimates$sate.aipw)


      if(bootstrap_zero){
        list(
          b = b,
          downwind_LogRain_b = downwind_LogRain_b,
          downwind_lmm_param_b = downwind_lmm_param_b,
          downwind_positive_target_lmm_param_b = downwind_positive_target_lmm_param_b,
          downwind_positive_control_lmm_param_b = downwind_positive_control_lmm_param_b,
          downwind_logistic_param_b = downwind_logistic_param_b,
          downwind_propensity_param_b = downwind_propensity_param_b,
          attr_b = attr_b,
          sate_b = sate_b
        )
      }else{
        list(
          b = b,
          downwind_LogRain_b = downwind_LogRain_b,
          downwind_lmm_param_b = downwind_lmm_param_b,
          downwind_positive_target_lmm_param_b = downwind_positive_target_lmm_param_b,
          downwind_positive_control_lmm_param_b = downwind_positive_control_lmm_param_b,
          attr_b = attr_b,
          sate_b = sate_b
        )
      }
    }

    parallel::stopCluster(cl)

    for(res in results){
      bootstrap_downwind_LogRain_matrix[res$b,] = res$downwind_LogRain_b
      bootstrap_downwind_lmm_param_matrix[res$b,] = res$downwind_lmm_param_b
      bootstrap_downwind_positive_target_lmm_param_matrix[res$b,] = res$downwind_positive_target_lmm_param_b
      bootstrap_downwind_positive_control_lmm_param_matrix[res$b,] = res$downwind_positive_control_lmm_param_b
      bootstrap_attr_matrix[res$b,] = res$attr_b
      bootstrap_sate_matrix[res$b,] = res$sate_b
      if(bootstrap_zero){
        bootstrap_downwind_logistic_param_matrix[res$b,] = res$downwind_logistic_param_b
        bootstrap_downwind_propensity_param_matrix[res$b,] = res$downwind_propensity_param_b
      }
    }

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

      bootstrap_downwind_propensity_param_matrix = bootstrap_downwind_propensity_param_matrix + matrix(data = rep( coef(ori_fitted_models$downwind_propensity_fit) - apply(bootstrap_downwind_propensity_param_matrix,2, function(x){mean(x, na.rm = T)}), B_bootstrap),
                                                                                                       nrow = B_bootstrap, ncol = ncol(bootstrap_downwind_propensity_param_matrix), byrow = TRUE)
    }



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
      downwind_LogRain = bootstrap_downwind_LogRain_matrix
    ))
  }else{
    return(list(
      hatattr = bootstrap_attr_matrix,
      hatsate = bootstrap_sate_matrix,
      downwind_lmm_param = bootstrap_downwind_lmm_param_matrix,
      downwind_logistic_param = NULL,
      downwind_propensity_param = NULL,
      downwind_positive_target_lmm_param = bootstrap_downwind_positive_target_lmm_param_matrix,
      downwind_positive_control_lmm_param = bootstrap_downwind_positive_control_lmm_param_matrix,
      downwind_LogRain = bootstrap_downwind_LogRain_matrix
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


#' @title Bootstrap Options
#' @description
#' This function generates a list of settings controlling how two-level bootstrap procedures are executed within \code{\link{rain_attr}}.
#'
#'
#' @param B_bootstrap An integer specifying the number of bootstrap replicates. Default is 10000.
#' @param bootstrap_type A character string specifying the type of bootstrap.
#'   Must be one of \code{"REB0"}, \code{"REB1"}, \code{"REB2"}, \code{"PREB0"},
#'   \code{"PREB1"}, \code{"PREB2"}, or \code{"MREB1"}. See \code{\link{bootstrap_downwind}} for their differences.
#'   Default is "PREB1".
#' @param bootstrap_zero Logical. If \code{TRUE}, the optional first-level bootstrap is performed to generate bootstrap samples of binary rainfall event indicators. Default is \code{TRUE}.
#' @param positive_prob_threshold An optional numeric value between 0 and 1 specifying the probability threshold for generating bootstrap samples of binary rainfall event indicators. Probabilities below this threshold are set to zero. Default is \code{NULL}.
#' @param discretize_rain Logical. If \code{TRUE}, rainfall values are discretized in bootstrap resamples. Default is \code{TRUE}.
#' @param winsorize_individual_rain Logical. If \code{TRUE}, individual rainfall values in bootstrap samples that exceed the upper bound specified by \code{individual_rain_interval} are replaced with random draws from a uniform distribution over \code{[individual_rain_interval[1], individual_rain_interval[2]]}. Default is \code{TRUE}.
#' @param individual_rain_interval Numeric vector of length 2 specifying the lower and upper bounds for adjusting bootstrapped individual rainfall values that are too large when \code{winsorize_individual_rain = TRUE}. Default is \code{c(100,175)}.
#' @param winsorize_total_rain Logical. If \code{TRUE}, all individual rainfall values in each bootstrap sample are proportionally rescaled so that the total equals a random number drawn uniformly from
#'   \code{[total_rain_interval[1], total_rain_interval[2]]} whenever the total bootstrapped rainfall falls outside this interval. Default is \code{TRUE}.
#' @param total_rain_interval Numeric vector of length 2 specifying the lower and upper bounds for adjusting the total of bootstrapped rainfall values when \code{winsorize_total_rain = TRUE}. Default is \code{c(6000,60000)}.
#' @param bootstrap_seed An integer specifying the random seed for the bootstrap procedure. Reproducibility is guaranteed only if \code{bootstrap_parallel} is the same, since parallel execution changes the order of random number generation. Default is 123.
#' @param bootstrap_parallel Logical. If \code{TRUE}, each bootstrap run is executed in parallel across multiple workers. If \code{FALSE}, they are run sequentially. Default is \code{FALSE}.
#' @param bootstrap_parallel_num_worker An integer specifying the number of parallel workers to use when \code{bootstrap_parallel = TRUE}. Default is \code{parallel::detectCores() - 1}.
#' @param CI_level A numeric value between 0 and 1 specifying the confidence level of the bootstrap percentile confidence intervals. Default is 0.95.
#'
#' @return A list containing all bootstrap options, suitable for passing to \code{\link{rain_attr}}.
#'
#' @details
#' This function is used to configure and store all settings needed for performing two-level bootstrap analyses within \code{\link{rain_attr}}.
#'
#' @seealso \code{\link{rain_attr}} for the main function, \code{\link{bootstrap_downwind}} for more details on the two-level bootstrap procedure
#'
#' @examples
#' #Create default bootstrap options to account for highly unbalanced clustered data
#' # Specifically: bootstrap_type = 'PREB1'
#' boot_options = bootstrap_option()
#' str(boot_options)
#'
#' #Bootstrap option with parallelization over (parallel::detectCores() - 1) number of workers and seed = 1 for reproducibility
#' boot_options_parallel = bootstrap_option(
#'   bootstrap_seed = 1,
#'   bootstrap_parallel = TRUE,
#'   bootstrap_parallel_num_worker = parallel::detectCores() - 1
#' )
#' str(boot_options_parallel)

bootstrap_option = function(B_bootstrap = 10000,
                            bootstrap_type = 'PREB1',
                            bootstrap_zero = T,
                            positive_prob_threshold = NULL,
                            discretize_rain = T,
                            winsorize_individual_rain = T,
                            individual_rain_interval = c(100,175),
                            winsorize_total_rain = T,
                            total_rain_interval = c(6000,60000),
                            bootstrap_seed = 123,
                            bootstrap_parallel = F,
                            bootstrap_parallel_num_worker = parallel::detectCores() - 1,
                            CI_level = 0.95){
  return(list(
    B_bootstrap = B_bootstrap,
    bootstrap_type = bootstrap_type,
    bootstrap_zero = bootstrap_zero,
    positive_prob_threshold = positive_prob_threshold,
    discretize_rain = discretize_rain,
    winsorize_individual_rain = winsorize_individual_rain,
    individual_rain_interval = individual_rain_interval,
    winsorize_total_rain = winsorize_total_rain,
    total_rain_interval = total_rain_interval,
    bootstrap_seed = bootstrap_seed,
    bootstrap_parallel = bootstrap_parallel,
    bootstrap_parallel_num_worker = bootstrap_parallel_num_worker,
    CI_level = CI_level
  ))
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
# attr_type = 'ThoEtAl'
# x_downwind_name = c('Gauge.Elevation', 'natural_pred')
# target_only = FALSE

#' @title Permutation-Based Procedure for Rainfall Enhancement Trial Data
#'
#' @description
#' Implements the permutation-based procedure used in \code{\link{rain_attr}} for inference on attribution and SATE in rainfall enhancement trial data under randomized assignment of ionizer operation.
#' It allows for flexible permutation schemes, including permutations between ionizers, between days, or between gauge-day combinations.
#'
#'
#' @param B_permutation An integer specifying the number of permutation replicates.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permute_between_ionizer Logical. If \code{TRUE}, for each day, a random permutation is performed among the operation statuses of all ionizers that have been deployed on that day.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permute_all_ionizers_between_day Logical. If \code{TRUE}, for each year, a random permutation is performed among the daily operation schedules of all trial days within that year.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permute_between_gaugeday Logical. If \code{TRUE}, for each year, a random permutation is performed among the gauge-day level operation schedule of all gauge-days within that year.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param ionizer_operation A data frame containing ionizer operation indicators for each day (row) and each ionizer (column), where 1 indicates that an ionizer is turned on and 0 indicates that it is off or not deployed yet.
#' Additionally, this data frame must include two columns with names specified by \code{ionizer_operation_day_column_name} and \code{ionizer_operation_year_column_name}, containing the day and year for each row.
#' Each day must appear only once in this data frame (no duplicated day entries).
#' The ionizer columns must appear in the same order as specified by \code{data_target_column_names} and must be consistent with the column order in \code{gaugeday_downwind}.
#'   (User-supplied using \code{\link{permutation_option}})
#' @param gaugeday_downwind A binary matrix indicating which gauge-day observations (row) are downwind of which ionizers (column), where 1 indicates that the gauge is downwind of the ionizer on that day, and 0 indicates that the gauge is not downwind of the ionizer or the ionizer has not been deployed yet.
#'   The row order of this matrix must match that of \code{data}. The column order must correspond to the ionizer columns in \code{ionizer_operation} (excluding the day and year columns) and be in the same order as specified by \code{data_target_column_names}.
#'   (User-supplied using \code{\link{permutation_option}})
#' @param year_ionizer_list A named list, where each element corresponds to a year and contains the names of deployed ionizers in that year.
#'   (User-supplied using \code{\link{permutation_option}})
#' @param data_target_column_names A character vector specifying the column names of \code{data} corresponding to the binary target indicators used in the downwind (second stage) LMM fitting.
#'   The order of names in this vector must match the column order of the corresponding ionizers in \code{ionizer_operation} (excluding the day and year columns) and in \code{gaugeday_downwind}.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param ionizer_operation_year_column_name A character string specifying the column name of \code{ionizer_operation} containing the year of each day.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param ionizer_operation_day_column_name A character string specifying the column name of \code{ionizer_operation} containing the day of each observation. The same column name should also be found in \code{data}.
#'   (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permutation_seed An integer specifying the random seed for the permutation-based procedure. Reproducibility is guaranteed only if \code{permutation_parallel} is the same, since parallel execution changes the order of random number generation.
#' (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permutation_parallel Logical. If \code{TRUE}, each permutation run is executed in parallel across multiple workers. If \code{FALSE}, they are run sequentially.
#' (User-configurable permutation option using \code{\link{permutation_option}})
#' @param permutation_parallel_num_worker An integer specifying the number of parallel workers to use when \code{permutation_parallel = TRUE}.
#' (User-configurable permutation option using \code{\link{permutation_option}})
#' @param data A data frame containing the original dataset used in \code{\link{rain_attr}}, along with an additional column containing the fitted values generated from the upwind (first stage) LMM.
#'   Its column names should contain \code{data_target_column_names} (binary target indicators) and \code{ionizer_operation_day_column_name} (day).
#'   The row order of this data frame must match that of \code{gaugeday_downwind}.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind_lmm_formula A two sided linear formula object to be used in \link[lme4]{lmer}, describing both the fixed-effects and random intercept part of the downwind (second stage) LMM.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind_propensity_formula A two sided linear formula object to be used in \code{\link{stats}{glm}} with \code{family = "binomial"}, for fitting a propensity score model to the treatment indicators of downwind (second stage) observations.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param attr_type A character string specifying the type of attribution estimates. Must be one of \code{"Chambers_Chandra"}, \code{"ThoEtAl"}, or \code{"No"}. See \code{\link{rain_attr}} for more information.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param x_downwind_name A character vector containing variable names from the right hand side of \code{downwind_lmm_formula}, for those variables that are not related to ionizers (treatment). The intercept is always included and does not need to be specified.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param target_only Logical. If \code{TRUE} the attribution estimates are computed based on only target observations. If \code{FALSE} the attribution estimates are computed based on both treatment and control observations.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param rain_col_name A character string specifying the column name of the raw scale rainfall in \code{ori_data}.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#'
#' @return A list with two components:
#' \describe{
#'   \item{hatattr}{Matrix of permutation replicates of attribution estimates.}
#'   \item{hatsate}{Matrix of permutation replicates of SATE estimates.}
#' }
#'
#' @details
#' This function implements the permutation-based procedure used in \code{\link{rain_attr}}. It is intended for internal use only. Users should not call this function directly.
#' Instead, permutation-based inference should be performed by calling \code{\link{rain_attr}} with \code{permutation = TRUE}.
#'
#  \strong{Additional Information Required for Permutation-Based Procedure} \cr
#' To perform permutation-based procedure, additional information need to be supplied through the following arguments of \code{\link{permutation_option}}:
#' \describe{
#'  \item{\code{year_ionizer_list}}{Timing of deployment of ionizers over the years.}
#'  \item{\code{ionizer_operation}}{Day(group)-level ionizers operation schedule during the rainfall enhancement trial.}
#'  \item{{gaugeday_downwind}}{Gauge-day(unit within group)-level information on relative orientation of gauges from ionizers each day.}
#' }
#'
#' \strong{Permutation Steps} \cr
#' The permutation-based procedure considers to randomly permute ionizers' operation statuses via:
#' \enumerate{
#'    \item{If \code{permute_between_ionizer = TRUE}, for each row of the day-level \code{ionizer_operation} matrix, a random permutation is performed among the binary indicators in the row that correspond to ionizers that have already been deployed during the year of the row.
#'    This is equivalent to randomly permuting the operation statuses of deployed ionizers for each day. }
#'    \item{If \code{permute_all_ionizers_between_day = TRUE}, for each year, a random permutation is performed among all rows belonging to the year in the day-level \code{ionizer_operation} matrix.
#'    This is equivalent to randomly permuting the daily operation schedules of all trial days belonging to each year.  }
#'    \item{The permuted day-level \code{ionizer_operation} are then expanded into a gauge-day level binary ionizers' operation indicator matrix, to match the gauge-day level \code{data}. }
#'    \item{If \code{permute_between_gaugeday = TRUE}, for each year, a random permutation is performed among all rows belonging to the year in the gauge-day level binary ionizers' operation indicator matrix.
#'    This is equivalent to randomly permuting the gauge-day level operation schedules of all gauge-days belonging to each year. }
#' }
#'
#' The above three optional permutation steps are performed in sequence, resulting in a final permuted gauge-day level binary ionizers' operation indicator matrix.
#' An elementwise multiplication is carried out between this matrix and \code{gaugeday_downwind}, and the results are used to replace the original columns (\code{data_target_column_names}) in \code{data} that contain the binary target indicators used in the downwind (second stage) LMM fitting.
#' Therefore, it is important to ensure that
#' \itemize{
#'  \item{The row order of \code{gaugeday_downwind} and \code{data} is consistent.}
#'  \item{The column orders of \code{gaugeday_downwind} and \code{ionizer_operation} is consistent, and matches the order specified by \code{data_target_column_names}.}
#' }
#' Based on these permuted binary target indicators, the original binary indicators \eqn{I_{ij}} for exposure to ionizers (treatment) are also updated to be their permuted counterparts \eqn{I_{ij}^*}, where \eqn{I_{ij}^*} only equals zero if all permuted binary target indicators in its corresponding row equal zero.
#'
#'
#'  \strong{Permutation Distribution of Attribution and SATE} \cr
#' The same procedure described in \code{\link{rain_attr}} is then repeated on the permuted dataset, where the columns containing the binary target indicators and the binary indicators \eqn{I_{ij}} for exposure to ionizers (treatment) are replaced according to the permutation.
#This includes the fitting of the downwind (second stage) LMM, downwind (second stage) target-only LMM, downwind (second stage) control-only LMM to the subset of observations from \code{data} that are downwind (second stage) and with positive rainfall, where a gauge-day level observation is said to be downwind if the gauge is downwind of at least one deployed ionizer (not necessarily turned on) on that day i.e., at least one of the indicators is equal to one for that row of \code{gaugeday_downwind}.
#' This includes the fitting of the downwind (second stage) LMM, downwind (second stage) target-only LMM, downwind (second stage) control-only LMM to the subset of observations from \code{data} that are downwind (second stage) and with positive rainfall, along with the fitting of downwind propensity score model to the subset of observations from \code{data} that are downwind (second stage) with the response being the permuted indicator \eqn{I_{ij}^*} for exposure to ionizers (treatment).
#' Finally, two attribution estimates and the SATE estimates are computed based on the estimation results of these models fitted to the permuted dataset, where \eqn{z_{ij}} (ionizer related covariate vector constructed from the binary target indicators) and \eqn{I_{ij}} are replaced by their permuted counterparts \eqn{z_{ij}^*} and \eqn{I_{ij}^*}, respectively.
#'
#' By repeatedly permuting ionizers' operation schedules, fitting models and computing attribution and SATE estimates for \code{B_permutation} number of times, this function returns the permutation distributions of
#' \itemize{
#'    \item{Two attribution estimates}
#'    \item{Five SATE estimates}
#' }
#'
#' \strong{Permutation-Based P-Value and Permutation-Based Plots} \cr
#' The permutation distributions of attribution and SATE produced by this function are further used in \code{\link{rain_attr}} to:
#' \itemize{
#'    \item{Compute permutation-based p-value as the proportion of permuted estimates that are greater than or equal to the original estimate, i.e.,
#'      \deqn{\frac{1}{B} \sum_{b=1}^{B} 1_{ \{ \hat{\theta}^*_{b} \geq \hat{\theta} \}  } ,}
#'      where \eqn{\hat{\theta}^*_{b}} is the estimate from the \eqn{b}-th permuted dataset and \eqn{\hat{\theta}} is the corresponding estimate from the original dataset.
#'    }
#'    \item{Plot the kernel density estimate with a solid vertical line for the original estimate \eqn{\hat{\theta}}. }
#' }
#'
#' @seealso \code{\link{rain_attr}} for the main function, \code{\link{permutation_option}} for specifying permutation options



permutation_ionizer = function(B_permutation, permute_between_ionizer, permute_all_ionizers_between_day, permute_between_gaugeday,
                               ionizer_operation, gaugeday_downwind, year_ionizer_list,
                               data_target_column_names, ionizer_operation_year_column_name, ionizer_operation_day_column_name,
                               permutation_seed, permutation_parallel, permutation_parallel_num_worker,
                               data, downwind_lmm_formula, downwind_propensity_formula,
                               attr_type, x_downwind_name, target_only,
                               rain_col_name){
  #data_target_column_names are the column names of 'data', which correspond to the target indicators of all ionizers
  #ionizer_operation_year_column_name is the column name of 'ionizer_operation' containing the year of each day
  #ionizer_operation_day_column_name is the column name of 'ionizer_operation' containing the day of each observation, which should be the same column name in 'data'
  #Note that data_target_column_names should have the same ordering as colnames(ionizer_operation), as well as gaugeday_downwind
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

  if(!permutation_parallel){
    set.seed(permutation_seed)

    downwind = apply(gaugeday_downwind,1,sum) > 0
    positive = ( data[,rain_col_name] > 0)

    test_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = data[downwind & positive,])
    z_downwind_name = setdiff(names(lme4::fixef(test_downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
    perm_downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)
    perm_downwind_propensity_formula = update.formula(downwind_propensity_formula, permuted_target_indicator ~ . )

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
        perm_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = perm_data[downwind & positive,])

        perm_downwind_positive_data = perm_data[downwind & positive,]
        target_vec = apply(perm_data[,data_target_column_names],1, sum) > 0
        nontarget_vec = apply(perm_data[,data_target_column_names],1, sum) ==  0
        perm_downwind_positive_target = target_vec[downwind & positive]
        perm_downwind_positive_control = nontarget_vec[downwind & positive]

        perm_hatattr = attr_est(attr_type, perm_downwind_positive_data, rain_col_name, perm_downwind_positive_target, perm_downwind_positive_control,
                                x_downwind_name, target_only = target_only, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)


        perm_downwind_positive_data$permuted_target_indicator = as.logical(perm_downwind_positive_target)


        perm_hatsate = sate_est(perm_downwind_positive_data, perm_downwind_positive_target, perm_downwind_positive_control, perm_downwind_propensity_formula, perm_downwind_separate_formula,
                                x_downwind_name, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

        perm_attr_matrix[b,] = c(perm_hatattr$apo, perm_hatattr$apl)
        perm_sate_matrix[b,] = c(perm_hatsate$estimates$sate.mb, perm_hatsate$estimates$sate.ipw, perm_hatsate$estimates$sate.ipw.l, perm_hatsate$estimates$sate.ipw.ma, perm_hatsate$estimates$sate.aipw)
      },error=function(e){cat(b,"th","Permutation Run Skipped due to ERROR :",conditionMessage(e), "\n")})
    }
  }else{
    downwind = apply(gaugeday_downwind,1,sum) > 0
    positive = ( data[,rain_col_name] > 0)

    test_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = data[downwind & positive,])
    z_downwind_name = setdiff(names(lme4::fixef(test_downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
    perm_downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)
    perm_downwind_propensity_formula = update.formula(downwind_propensity_formula, permuted_target_indicator ~ . )

    permutation_cl = parallel::makeCluster(permutation_parallel_num_worker)
    parallel::clusterExport(permutation_cl, varlist = c('attr_est', 'sate_est'))
    doParallel::registerDoParallel(permutation_cl)
    doRNG::registerDoRNG(seed = permutation_seed)
    `%dopar%` <- foreach::`%dopar%`


    permutation_results = foreach::foreach(b = 1:B_permutation, .packages = c("dplyr","lme4"), .errorhandling = 'remove') %dopar% {
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
      perm_downwind_lmm_fit = lme4::lmer(downwind_lmm_formula, data = perm_data[downwind & positive,])

      perm_downwind_positive_data = perm_data[downwind & positive,]
      target_vec = apply(perm_data[,data_target_column_names],1, sum) > 0
      nontarget_vec = apply(perm_data[,data_target_column_names],1, sum) ==  0
      perm_downwind_positive_target = target_vec[downwind & positive]
      perm_downwind_positive_control = nontarget_vec[downwind & positive]

      perm_hatattr = attr_est(attr_type, perm_downwind_positive_data, rain_col_name, perm_downwind_positive_target, perm_downwind_positive_control,
                              x_downwind_name, target_only = target_only, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)


      perm_downwind_positive_data$permuted_target_indicator = as.logical(perm_downwind_positive_target)



      perm_hatsate = sate_est(perm_downwind_positive_data, perm_downwind_positive_target, perm_downwind_positive_control, perm_downwind_propensity_formula, perm_downwind_separate_formula,
                              x_downwind_name, downwind_lmm_fit = perm_downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)

      perm_attr_b = c(perm_hatattr$apo, perm_hatattr$apl)
      perm_sate_b = c(perm_hatsate$estimates$sate.mb, perm_hatsate$estimates$sate.ipw, perm_hatsate$estimates$sate.ipw.l, perm_hatsate$estimates$sate.ipw.ma, perm_hatsate$estimates$sate.aipw)

      list(
        b = b,
        perm_attr_b = perm_attr_b,
        perm_sate_b = perm_sate_b
      )
    }

    parallel::stopCluster(permutation_cl)

    for(permutation_res in permutation_results){
      perm_attr_matrix[permutation_res$b,] = permutation_res$perm_attr_b
      perm_sate_matrix[permutation_res$b,] = permutation_res$perm_sate_b
    }

  }


  return(list(
    hatattr = perm_attr_matrix,
    hatsate = perm_sate_matrix
  ))
}


#' @title Permutation Options
#'
#' @description
#' This function generates a list of settings controlling how permutation-based procedures are executed within \code{\link{rain_attr}}.
#'
#'
#' @param B_permutation An integer specifying the number of permutation replicates. Default is 10000.
#' @param permute_between_ionizer Logical. If \code{TRUE}, for each day, a random permutation is performed among the operation statuses of all ionizers that have been deployed on that day. Default is \code{TRUE}.
#' @param permute_all_ionizers_between_day Logical. If \code{TRUE}, for each year, a random permutation is performed among the daily operation schedules of all trial days within that year. Default is \code{FALSE}.
#' @param permute_between_gaugeday Logical. If \code{TRUE}, for each year, a random permutation is performed among the gauge-day level operation schedule of all gauge-days within that year. Default is \code{TRUE}.
#' @param ionizer_operation_input A data frame containing ionizer operation indicators for each day (row) and each ionizer (column), where 1 indicates that an ionizer is turned on and 0 indicates that it is off or not deployed yet.
#' Additionally, this data frame must include two columns with names specified by \code{ionizer_operation_day_column_name} and \code{ionizer_operation_year_column_name}, containing the day and year for each row.
#' Each day must appear only once in this data frame (no duplicated day entries).
#' The ionizer columns must appear in the same order as specified by \code{data_target_column_names} and must be consistent with the column order in \code{gaugeday_downwind_input}. Default is \code{ionizer_operation}.
#' @param gaugeday_downwind_input A binary matrix indicating which gauge-day observations (row) are downwind of which ionizers (column), where 1 indicates that the gauge is downwind of the ionizer on that day, and 0 indicates that the gauge is not downwind of the ionizer or the ionizer has not been deployed yet.
#'   The row order of this matrix must match that of the original dataset supplied to \code{\link{rain_attr}}. The column order must correspond to the ionizer columns in \code{ionizer_operation_input} (excluding the day and year columns) and be in the same order as specified by \code{data_target_column_names}. Default is \code{gaugeday_downwind}.
#' @param year_ionizer_list A named list, where each element corresponds to a year and contains the names of deployed ionizers in that year. Default is:
#'   \code{list('2013' = c('H1','H2'), '2014' = c('H1','H2','H3','H4'), '2015' = c('H1','H2','H3','H4','H5','H6'), '2016' = c('H1','H2','H3','H4','H5','H6','H7','H8'), '2017' = c('H1','H2','H3','H4','H5','H6','H7','H8','H9','H10'), '2018' = c('H1','H2','H3','H4','H5','H6','H7','H8','H9','H10'))}.
#' @param data_target_column_names A character vector specifying the column names of the original dataset supplied to \code{\link{rain_attr}}, corresponding to the binary target indicators used in the downwind (second stage) LMM fitting.
#'   The order of names in this vector must match the column order of the corresponding ionizers in \code{ionizer_operation_input} (excluding the day and year columns) and in \code{gaugeday_downwind_input}. Default is:
#'   \code{c("Target.H.01", "Target.H.02", "Target.H.03", "Target.H.04", "Target.H.05", "Target.H.06", "Target.H.07", "Target.H.08", "Target.H.09", "Target.H.10")}.
#' @param ionizer_operation_year_column_name A character string specifying the column name of \code{ionizer_operation_input} containing the year of each day. Default is \code{'Year'}.
#' @param ionizer_operation_day_column_name A character string specifying the column name of \code{ionizer_operation_input} containing the day of each observation. The same column name should also be found in the original dataset supplied to \code{\link{rain_attr}}. Default is \code{'TrialDay'}.
#' @param permutation_seed An integer specifying the random seed for the permutation-based procedure. Reproducibility is guaranteed only if \code{permutation_parallel} is the same, since parallel execution changes the order of random number generation. Default is 123.
#' @param permutation_parallel Logical. If \code{TRUE}, each permutation run is executed in parallel across multiple workers. If \code{FALSE}, they are run sequentially. Default is \code{FALSE}.
#' @param permutation_parallel_num_worker An integer specifying the number of parallel workers to use when \code{permutation_parallel = TRUE}. Default is \code{parallel::detectCores() - 1}.
#'
#' @return A list containing all permutation options, suitable for passing to \code{\link{rain_attr}}.
#'
#' @details
#' This function is used to configure and store all settings needed for performing permutation-based analyses within \code{\link{rain_attr}}.
#' It is important to ensure that
#' \itemize{
#'  \item{The row order of \code{gaugeday_downwind_input} and the original dataset supplied to \code{\link{rain_attr}} is consistent.}
#'  \item{The column orders of \code{gaugeday_downwind_input} and \code{ionizer_operation_input} is consistent, and matches the order specified by \code{data_target_column_names}.}
#' }
#'
#' @seealso \code{\link{rain_attr}} for the main function, \code{\link{permutation_ionizer}} for more details on the permutation-based procedure
#'
#' @examples
#' #Create default permutation options
#' # These are the same permutation settings used in Chambers et al. (2022)
#' # "Nudging a Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall Enhancement Trial in Oman"
#' # Specifically: permute_between_ionizer = TRUE, permute_all_ionizers_between_day = FALSE, and permute_between_gaugeday = TRUE
#' perm_options = permutation_option()
#' str(perm_options)
#'
#' #Permutation option with parallelization over (parallel::detectCores() - 1) number of workers and seed = 1 for reproducibility
#' perm_options_parallel = permutation_option(
#'   permutation_seed = 1,
#'   permutation_parallel = TRUE,
#'   permutation_parallel_num_worker = parallel::detectCores() - 1
#' )
#' str(perm_options_parallel)
#'
permutation_option = function(B_permutation = 10000,
                              permute_between_ionizer = T,
                              permute_all_ionizers_between_day = F,
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
                              permutation_parallel_num_worker = parallel::detectCores() - 1){
  return(list(
    B_permutation = B_permutation,
    permute_between_ionizer = permute_between_ionizer,
    permute_all_ionizers_between_day = permute_all_ionizers_between_day,
    permute_between_gaugeday = permute_between_gaugeday,
    ionizer_operation = ionizer_operation_input,
    gaugeday_downwind = gaugeday_downwind_input,
    year_ionizer_list = year_ionizer_list,
    data_target_column_names = data_target_column_names,
    ionizer_operation_year_column_name = ionizer_operation_year_column_name,
    ionizer_operation_day_column_name = ionizer_operation_day_column_name,
    permutation_seed = permutation_seed,
    permutation_parallel = permutation_parallel,
    permutation_parallel_num_worker = permutation_parallel_num_worker
  ))
}

permutation_p_value = function(permutation_result, ori_est){
  temp = sapply(1:ncol(permutation_result), FUN = function(i){
    mean(permutation_result[,i] >= ori_est[i], na.rm = T )
  })
  names(temp) = colnames(permutation_result)
  return(temp)
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




#For checking: apo should be 0.111206, apl should be 0.1251201 for attr_type = 'Chambers_Chandra'
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
#                   bootstrap_option = bootstrap_option(B_bootstrap = 3,
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
#                         bootstrap_option = bootstrap_option(B_bootstrap = 3,
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
#                         bootstrap_option = bootstrap_option(B_bootstrap = 10,
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
#                        bootstrap_option = bootstrap_option(B_bootstrap = 10,
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
#                                       bootstrap_option = bootstrap_option(B_bootstrap = 10,
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
#                                                bootstrap_option = bootstrap_option(B_bootstrap = 10,
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
#                     permutation_option = permutation_option(
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
# my_perm_result_TT_Chambers_Chandra = rain_attr(data = oman,
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
#                                            attr_type = 'Chambers_Chandra',
#                                            x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                            target_only = FALSE,
#                                            bootstrap =F,
#                                            bootstrap_option = NULL,
#                                            permutation = T,
#                                            permutation_option = permutation_option(
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
# max(abs(perm_result_TT$perm_attribution_Chambers_Chandra_matrix[1:6,c('apo','apl')] - my_perm_result_TT_Chambers_Chandra$permutation_result$hatattr))
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
#                                            permutation_option = permutation_option(
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
# my_perm_result_TT_Chambers_Chandra = rain_attr(data = oman,
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
#                                            attr_type = 'Chambers_Chandra',
#                                            x_downwind_name = c('Gauge.Elevation', 'natural_pred'),
#                                            target_only = FALSE,
#                                            bootstrap =F,
#                                            bootstrap_option = NULL,
#                                            permutation = T,
#                                            permutation_option = permutation_option(
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
# max(abs(perm_result_TT$perm_attribution_Chambers_Chandra_matrix[1:6,c('apo','apl')] - my_perm_result_TT_Chambers_Chandra$permutation_result$hatattr))
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
#                                        permutation_option = permutation_option(
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
