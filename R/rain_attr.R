#' @title
#' Attribution and Sample Average Treatment Effect for Rainfall Enhancement Trial Data
#'
#' @description
#' Perform estimation and inference of attribution and sample average treatment effect for rainfall enhancement trial data, based on the two-stage linear mixed model (LMM) approach employed in Chambers et al. (2022a).
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
#' \item{\code{ChambersEtAl}}{Attribution is estimated based on the approach of Chambers et al. (2022a), to adjust for back-transformation bias due to the modelling of log-transformed rainfall:
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
#' \item{\code{ChambersEtAl_No_Winsorize}}{Similar to \code{ChambersEtAl}, but without the winsorizing step, i.e., \eqn{\max\{\lambda^{-1} \exp(-z_{ij}^\top \hat{\beta}), 0.5\}} is replaced by \eqn{\lambda^{-1} \exp(-z_{ij}^\top \hat{\beta})}:
#'     \deqn{
#'     \code{apo} = \sum_{(i,j)} Rain_{ij} \{ 1 - \lambda^{-1} \exp(-z_{ij}^\top \hat{\beta}) \} /  \sum_{(i,j)}Rain_{ij}, \quad
#'     \code{apl} = \sum_{(i,j)} Rain_{ij} \{ 1 - \lambda^{-1} \exp(-z_{ij}^\top \hat{\beta}) \} /  \sum_{(i,j)}Rain_{ij} \lambda^{-1} \exp(-z_{ij}^\top \hat{\beta}),
#'     }
#'     where \eqn{\lambda} is the same as in \code{ChambersEtAl}.
#'   }
#'
#' \item{\code{ThoEtAl}}{Attribution is estimated based on an alternative adjustment proposed by Tho et al. (2026), using the estimated covariance matrix \eqn{\hat{\Sigma}} of \eqn{\hat{\beta}}.
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
#' These IPW weights are then used, together with the estimation results of the downwind (second stage) LMM, to obtain the following five types of SATE estimates discussed in Chambers et al. (2022b):
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
#' This function can also be used to perform bootstrap inference on the attribution and SATE, by setting \code{bootstrap = TRUE} and supplying the relevant bootstrap options using \code{\link{bootstrap_opt}()}.
#' For full details of the bootstrap procedure, please see \code{\link{bootstrap_downwind}}. Briefly, the bootstrap is carried out in two levels by conditioning on the upwind (first stage) LMM and its fitted values:
#' \itemize{
#' \item{First level is an optional level that is only carried out when \code{bootstrap_opt(bootstrap_zero = TRUE)}. This level considers generating bootstrap samples of rainfall event indicator for the subset of observations satisfying \code{downwind_subset} using the predicted probabilities from the downwind logistic model. This model is fitted using \code{glm(downwind_logistic_formula, family = "binomial")} to the subset of observations from \code{data} satisfying \code{downwind_subset}, with the response being an indicator for whether the observed rainfall is greater than zero, i.e., the indicator is defined by the logical expression \code{positive_subset}.
#'   }
#'
#' \item{Second level generates bootstrap samples of positive rainfall for the subset of observations not only satisfying \code{downwind_subset} but also with the first-level bootstrapped rainfall event indicator being equal to one. When \code{bootstrap_opt(bootstrap_zero = FALSE)}, then this level generates bootstrap samples of positive rainfall for the subset of observations satisfying \code{downwind_subset & positive_subset}.
#'   This is done using one of the semiparametric bootstrap methods of Chambers & Chandra (2013) and Tho et al. (2025), which involves the use of marginal residuals from the fitted downwind (second stage) LMM.   }
#' }
#' The above attribution and SATE estimates are then computed based on each bootstrap sample of the positive rainfall, forming their respective bootstrap distributions. This function also provides bootstrap distributions of parameters associated with the downwind LMM (\code{downwind_lmm_formula}), downwind logistic model (\code{downwind_logistic_formula}), downwind propensity score model (\code{downwind_propensity_formula}), downwind treatment-only LMM, and downwind control-only LMM.
#' These bootstrap distributions are then used to compute bootstrap p-values (proportion of bootstrapped estimates that are negative), form bootstrap percentile confidence intervals (with confidence level specified in \code{bootstrap_option$CI_level}), and generate their respective plots.
#' It is worth noting that the entire bootstrap procedure (including rainfall resampling, model fitting, and parameter estimation) can be run in parallel by setting \code{bootstrap_option$bootstrap_parallel = TRUE}, using \code{bootstrap_option$bootstrap_parallel_num_worker} workers.
#'
#' Finally, this function enables permutation-based inference on the attribution and SATE, by setting \code{permutation = TRUE} and supplying the relevant permutation options using \code{\link{permutation_opt}()}.
#' For full details of the permutation-based procedure, please see \code{\link{permutation_ionizer}}. In short, the permutation-based procedure involves randomly permuting the operating schedules of the ionizers (treatment) and re-estimating the attribution and SATE based on the permuted data, from which permutations distributions of attribution and SATE estimates are formed.
#' These permutation distributions are used to compute permutation p-values (proportion of permuted estimates that are greater than the observed estimates) and generate their respective plots.
#'
#'
#' @param data A data frame containing the variables named in \code{upwind_lmm_formula}, \code{downwind_lmm_formula}, \code{downwind_logistic_formula} (if specified), and \code{downwind_propensity_formula}.
#' It should also contain variables named in \code{rain_col_name}, \code{upwind_subset}, \code{downwind_subset}, \code{downwind_target_subset}, and \code{downwind_control_subset}.
#' @param upwind_lmm_formula A two sided linear formula object to be used in \link[lme4]{lmer}, describing both the fixed-effects and random intercept part of the upwind (first stage) LMM.
#' @param instr_pred_name A character string to store the variable name of the fitted values generated from the upwind (first stage) LMM.
#' @param instr_pred_type Type of fitted values generated from the upwind (first stage) LMM. If "Unconditional" the fitted values equal to only the estimated fixed effects. If "Conditional" the fitted values equal to the sum of estimated fixed effects and EBLUPs of random intercepts.
#' @param downwind_lmm_formula A two sided linear formula object to be used in \link[lme4]{lmer}, describing both the fixed-effects and random intercept part of the downwind (second stage) LMM. This formula should contain the variable name specified in \code{instr_pred_name}.
#' @param downwind_logistic_formula An optional two sided linear formula object to be used in \code{\link{glm}} with \code{family = "binomial"}, for fitting a logistic model to the indicators of rainfall event. This only needs to be specified when \code{bootstrap = TRUE} and \code{bootstrap_option$bootstrap_zero = TRUE}.
#' @param downwind_propensity_formula A two sided linear formula object to be used in \code{\link{glm}} with \code{family = "binomial"}, for fitting a propensity score model to the treatment indicators of downwind (second stage) observations.
#' @param rain_col_name A character string that refers to the column name of the raw scale rainfall in \code{data}.
#' @param upwind_subset A logical expression used to extract the relevant subset of observations from \code{data} to be used in the upwind (first stage) LMM fitting. For example, \code{Gauge.Day.Type == "Upwind"}.
#' @param downwind_subset A logical expression used to extract the relevant subset of observations from \code{data} to be used in the downwind (second stage) LMM fitting. For example, \code{Gauge.Day.Type \%in\% c("Target","Control")}.
#' @param downwind_target_subset A logical expression used to extract the relevant subset of downwind (second stage) observations from \code{data} that were exposed to treatment (operating ionizers). For example, \code{Gauge.Day.Type == "Target"}.
#' @param downwind_control_subset A logical expression used to extract the relevant subset of downwind (second stage) observations from \code{data} that were not exposed to treatment (operating ionizers). For example, \code{Gauge.Day.Type == "Control"}.
#' @param positive_subset A logical expression used to extract the relevant subset of observations from \code{data} with positive rainfall - these are the observations that are used in the fitting of upwind (first stage) LMM, downwind (second stage) LMM, downwind (second stage) treatment-only LMM, downwind (second stage) control-only LMM, and the downwind (second stage) propensity score model.
#' @param attr_type An optional character string specifying the type of attribution estimates. Must be one of \code{"ChambersEtAl"}, \code{"ChambersEtAl_No_Winsorize"}, \code{"ThoEtAl"} (default), or \code{"No"}. See "Details" for more information.
#' @param x_downwind_name A character vector containing variable names from the right hand side of \code{downwind_lmm_formula}, for those variables that are not related to ionizers (treatment). The intercept is always included and does not need to be specified.
#' @param target_only An optional logical. If \code{TRUE} the attribution estimates are computed based on only treated observations. If \code{FALSE} the attribution estimates are computed based on both treated and control observations.
#' @param bootstrap An optional logical. If \code{TRUE} bootstrap is carried out to perform inference on the attribution and sample average treatment effect. If \code{FALSE} (default) no bootstrap is carried out.
#' @param bootstrap_option An optional list containing all bootstrap settings, used only when \code{bootstrap = TRUE}. See \code{\link{bootstrap_opt}} for the default list elements and their usage.
#' @param permutation An optional logical, If \code{TRUE} randomized permutation is carried out on the ionizer operation (treatment) schedule to perform inference on the attribution and sample average treatment effect. If \code{FALSE} (default) no randomized permutation is carried out.
#' @param permutation_option An optional list containing all permutation settings, used only when \code{permutation = TRUE}. See \code{\link{permutation_opt}} for the default list elements and their usage.
#'
#'
#' @returns An object of class \code{\link[=rain_attr-class]{rain_attr}}, which is a list containing
#' \describe{
#' \item{all_fitted_models}{A list of model objects from \code{\link[lme4]{lmer}} for the first stage (upwind), second stage (downwind) LMM, second stage (downwind) treatment-only LMM, and second stage (downwind) control-only LMM, along with model objects from \code{\link{glm}} for the logistic model of rainfall event indicator (\code{NULL} if \code{downwind_logistic_formula} is not specified) and the propensity score model for the treatment indicator of second stage (downwind) observations.}
#' \item{hatattr}{A vector containing the attribution estimates.}
#' \item{hatsate}{A vector containing the sample average treatment effect estimates.}
#' \item{bootstrap_result}{A list of matrices with the following elements:
#'
#'  - hatattr: Matrix of bootstrap samples for attribution estimates.
#'  - hatsate: Matrix of bootstrap samples for SATE estimates.
#'  - downwind_lmm_param: Matrix of bootstrap samples for fixed effect coefficient and random effect variance estimates of downwind (second stage) LMM.
#'  - downwind_logistic_param: Matrix of bootstrap samples for regression coefficient estimates of downwind logistic model fitted to the rainfall event indicators. This is \code{NULL} if \code{downwind_logistic_formula} is not specified.
#'  - downwind_propensity_param: Matrix of bootstrap samples for regression coefficient estimates of downwind propensity score model fitted to the treatment indicators.
#'  - downwind_positive_target_lmm_param: Matrix of bootstrap samples for fixed effect coefficient and random effect variance estimates of downwind (second stage) treatment-only LMM.
#'  - downwind_positive_control_lmm_param: Matrix of bootstrap samples for fixed effect coefficient and random effect variance estimates of downwind (second stage) control-only LMM.
#'  - downwind_LogRain: Matrix of bootstrap samples for the log-transformed rainfall of all downwind (second-stage) observations. Observations with zero bootstrapped rainfall are represented as \code{NA}.
#' }
#'
#' \item{bootstrap_CI_result}{A list of matrices with same element names as in \code{bootstrap_result} (excluding \code{downwind_LogRain}), containing the corresponding bootstrap percentile confidence intervals.}
#' \item{bootstrap_p_value_result}{A list of matrices with same element names as in \code{bootstrap_result} (excluding \code{downwind_LogRain}), containing the corresponding proportion of bootstrap samples that are less than zero.}
#' \item{bootstrap_plot_result}{A list of matrices with two elements:
#'
#' - hatattr: A list of \code{ggplot} objects, each showing the bootstrap distribution of attribution estimates. Each plot includes a dotted vertical line at zero and a solid vertical line at the original estimate based on the observed data.
#' - hatsate: A list of \code{ggplot} objects, each showing the bootstrap distribution of SATE estimates. Each plot includes a dotted vertical line at zero and a solid vertical line at the original estimate based on the observed data.
#' }
#'
#' \item{permutation_result}{A list of matrices with the two elements:
#'
#' - hatattr: Matrix of permutation samples for attribution estimates.
#' - hatsate: Matrix of permutation samples for SATE estimates.
#' }
#'
#' \item{permutation_p_value_result}{A list of matrices with same element names as in \code{permutation_result}, containing the corresponding proportion of permutation samples that are greater than or equal to the original estimate based on the observed data.}
#' \item{permutation_plot_result}{A list of matrices with two elements:
#'
#' - hatattr: A list of \code{ggplot} objects, each showing the permutation distribution of attribution estimates. Each plot includes a solid vertical line at the original estimate based on the observed data.
#' - hatsate: A list of \code{ggplot} objects, each showing the permutation distribution of SATE estimates. Each plot includes a solid vertical line at the original estimate based on the observed data.
#' }
#'
#' \item{args}{A list of the original function arguments.}
#' \item{data}{A data frame containing the original supplied \code{data}, with an additional column containing the fitted values generated from the upwind (first stage) LMM.}
#'
#'}
#'
#'@references
#'\itemize{
#'  \item Chambers, R., Beare, S., Peak, S. and Al-Kalbani, M. (2022a) Nudging a Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall Enhancement Trial in Oman. \emph{International Statistical Review}, 90: 346–373.
#'  \item Chambers, R. and Chandra, H. (2013). A Random Effect Block Bootstrap for Clustered Data. \emph{Journal of Computational and Graphical Statistics}, 22, 452–470.
#'  \item Chambers, R., Ranjbar, S., Salvati, N., and Pacini, B. (2022b) Weighting, Informativeness and Causal Inference, with an Application to Rainfall Enhancement. \emph{Journal of the Royal Statistical Society Series A: Statistics in Society}, 185: 1584–1612
#'  \item Tho, Z. Y., Chambers, R., and Welsh, A. H. (2025) A Proportional Random Effect Block Bootstrap for General Clustered Data.
#'     \href{https://arxiv.org/abs/2510.07770}{arXiv:2510.07770}.
#'  \item Tho, Z. Y., Chambers, R., and Welsh, A. H. (2026) Bias-Adjusted Attribution Estimation for Rainfall Enhancement Trials.
#'}
#'
#'@export

rain_attr = function(data, upwind_lmm_formula, instr_pred_name, instr_pred_type,
                     downwind_lmm_formula, downwind_logistic_formula = NULL, downwind_propensity_formula,
                     rain_col_name,
                     upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset, positive_subset,
                     attr_type = 'ThoEtAl', x_downwind_name, target_only = FALSE,
                     bootstrap = FALSE, bootstrap_option = bootstrap_opt(),
                     permutation = FALSE, permutation_option = permutation_opt()
                     ){

  if(!instr_pred_name %in% all.vars(downwind_lmm_formula)){
    stop("instr_pred_name cannot be found in downwind_lmm_formula")
  }

  if(formula.tools::lhs(stats::as.formula(gsub("[()]", "", downwind_propensity_formula))) != substitute(downwind_target_subset)){
    stop("The definition of Target provided in  downwind_target_subset is not consistent with the LHS of downwind_propensity_formula")
  }

  if(!is.null(downwind_logistic_formula)){
    if(formula.tools::lhs(stats::as.formula(gsub("[()]", "", downwind_logistic_formula))) != substitute(positive_subset)){
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

  if(!attr_type %in% c('ChambersEtAl', 'ChambersEtAl_No_Winsorize','ThoEtAl', 'No')){
    stop("attr_type should be one of 'ChambersEtAl', 'ChambersEtAl_No_Winsorize,'ThoEtAl', or 'No'")
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
    if(!permutation_option$ionizer_operation_year_column_name %in% colnames(permutation_option$ionizer_operation_input) ){
      stop("The columns of permutation_option$ionizer_operation_input do not contain permutation_option$ionizer_operation_year_column_name")
    }

    if(!permutation_option$ionizer_operation_day_column_name %in% colnames(permutation_option$ionizer_operation_input) ){
      stop("The columns of permutation_option$ionizer_operation_input do not contain permutation_option$ionizer_operation_day_column_name")
    }

    if(!permutation_option$ionizer_operation_day_column_name %in% colnames(data) ){
      stop("permutation_option$ionizer_operation_day_column_name cannot be found in the column names of data")
    }

    if(mean(permutation_option$data_target_column_names %in% colnames(data)) != 1){
      stop("At least one element of permutation_option$data_target_column_names cannot be found in the column names of data")
    }

    if(length(unique(permutation_option$ionizer_operation_input[,permutation_option$ionizer_operation_day_column_name])) != length(permutation_option$ionizer_operation_input[,permutation_option$ionizer_operation_day_column_name]) ){
      stop('permutation_option$ionizer_operation_input[,permutation_option$ionizer_operation_day_column_name] contains duplicated values')
    }

    if(mean(data[,permutation_option$ionizer_operation_day_column_name] %in% permutation_option$ionizer_operation_input[,permutation_option$ionizer_operation_day_column_name] ) !=1 ){
      stop('At least one day in data cannot be found in permutation_option$ionizer_operation_input')
    }

    if(nrow(data)!= nrow(permutation_option$gaugeday_downwind_input)){
      stop('Number of rows in data does not equal number of rows in permutation_option$gaugeday_downwind_input')
    }


  }

  original_args <- as.list(match.call())[-1]

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
  #Compute Point Estimates for Attribution - using ChambersEtAl or ChambersEtAl_No_Winsorize or ThoEtAl or No Estimates
  hatattr = attr_est(attr_type, downwind_positive_data, rain_col_name, downwind_positive_target, downwind_positive_control,
                     x_downwind_name, target_only = target_only, downwind_lmm_fit = fitted_models$downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)




  #Compute Points Estimates for SATE - using different types of SATE estimates
  z_downwind_name = setdiff(names(lme4::fixef(fitted_models$downwind_lmm_fit)), c('(Intercept)',x_downwind_name))
  downwind_separate_formula = remove_fixed_terms(input_formula = downwind_lmm_formula, vars_to_remove = z_downwind_name)
  hatsate = sate_est(downwind_positive_data, downwind_positive_target, downwind_positive_control, downwind_propensity_formula, downwind_separate_formula,
                     x_downwind_name, downwind_lmm_fit = fitted_models$downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL)




  all_fitted_models = list(
    upwind_lmm_fit = fitted_models$upwind_lmm_fit,
    downwind_lmm_fit = fitted_models$downwind_lmm_fit,
    downwind_logistic_fit = fitted_models$downwind_logistic_fit,
    downwind_propensity_fit = hatsate$fitted_models$downwind_propensity_fit,
    downwind_positive_target_lmm_fit = hatsate$fitted_models$downwind_positive_target_lmm_fit,
    downwind_positive_control_lmm_fit = hatsate$fitted_models$downwind_positive_control_lmm_fit
  )

  #Perform Bootstrap Inference

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




  #Perform Permutation Inference

  if(permutation){
    permutation_result = permutation_ionizer(B_permutation = permutation_option$B_permutation,
                                             permute_between_ionizer = permutation_option$permute_between_ionizer,
                                             permute_all_ionizers_between_day = permutation_option$permute_all_ionizers_between_day,
                                             permute_between_gaugeday = permutation_option$permute_between_gaugeday,
                                             ionizer_operation_input = permutation_option$ionizer_operation_input,
                                             gaugeday_downwind_input = permutation_option$gaugeday_downwind_input,
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

  output = list(
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
    args = original_args,
    data = fitted_models$data
  )

  class(output) = "rain_attr"

  return(output)
}




attr_est = function(attr_type, downwind_positive_data, rain_col_name, downwind_positive_target, downwind_positive_control,
                    x_downwind_name, target_only, downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL){
  #Note that hatu should be for downwind positive (i,t) regardless of target_only
  #Optional input arguments: hatalphabeta, hatu - these are only placeholders for now, which would be useful in the future if we want to develop alternative variant of REB2/PREB2 that considers mean-correcting variance components and hatalphabeta, before plugging them into the attribution estimate formula
  #Note that for attr_type == 'ThoEtAl', the hatSigma_beta matrix is ALWAYS obtained from downwind_lmm_fit
  #When we consider hatbeta from MQ, can add another option to use either hatbeta_{0.5} or hatbeta_{hat{q}_i}

  if(target_only){
    downwind_positive_useful_row = downwind_positive_target
  }else{
    downwind_positive_useful_row = (downwind_positive_target | downwind_positive_control)
  }

  y_vec = downwind_positive_data[downwind_positive_useful_row,rain_col_name]
  x_z_mat = stats::model.matrix(downwind_lmm_fit, data = downwind_positive_data, type = 'fixed')[downwind_positive_useful_row,]

  if(is.null(hatalphabeta)){
    hatalpha_downwind = lme4::fixef(downwind_lmm_fit)[c('(Intercept)',x_downwind_name)]
    hatbeta_downwind = lme4::fixef(downwind_lmm_fit)[setdiff(names(lme4::fixef(downwind_lmm_fit)), c('(Intercept)',x_downwind_name))]
  }else{
    hatalpha_downwind = hatalphabeta[c('(Intercept)',x_downwind_name)]
    hatbeta_downwind = hatalphabeta[setdiff(names(hatalphabeta), c('(Intercept)',x_downwind_name))]
  }



  if(attr_type == 'ChambersEtAl'){
    if(is.null(hatu)){
      hatu = stats::predict(downwind_lmm_fit, newdata = downwind_positive_data, random.only = TRUE)
    }

    #compute log_hatw = (1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t, for PDR (i,t) or PDR Target (i,t)
    #Depends on if there is offset term on LHS of formula(downwind_lmm_fit)
    if(length(formula.tools::lhs.vars(stats::formula(downwind_lmm_fit))) == 1){
      log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row])
    }

    if(length(formula.tools::lhs.vars(stats::formula(downwind_lmm_fit))) > 1){
      all_offset_terms = formula.tools::lhs.vars(stats::formula(downwind_lmm_fit))[2:length(formula.tools::lhs.vars(stats::formula(downwind_lmm_fit)))]
      if(length(all_offset_terms) > 1){
        log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row]) + apply(downwind_positive_data[,all_offset_terms],1,sum)
      }

      if(length(all_offset_terms) == 1){
        log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row]) + as.vector(downwind_positive_data[,all_offset_terms])
      }
    }


    #compute log_hatd = z_it %*% hatbeta, for PDR (i,t) or PDR Target (i,t)
    log_hatd = as.vector(x_z_mat[,setdiff(colnames(x_z_mat),c('(Intercept)',x_downwind_name))] %*% hatbeta_downwind )

    #compute log_haty_naive = (1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t + z_it %*% hatbeta, for PDR (i,t) or PDR Target (i,t)
    log_haty_naive = log_hatw + log_hatd

    #compute mu = (1/n) * sum of  (y_{i't'} / haty_{i't'}^naive ), where the sum is across PDR (i,t) or PDR Target (i,t), and n is the number of PDR (i,t) or PDR Target (i,t)  observation.
    mu = mean( y_vec / exp(log_haty_naive) )


    #compute m = Var{(1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t} / Var{z_it %*% hatbeta}, where the Var is across PDR (i,t) or PDR Target (i,t).
    m = (stats::var(log_hatw))/stats::var(log_hatd)

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

  if(attr_type == 'ChambersEtAl_No_Winsorize'){
    if(is.null(hatu)){
      hatu = stats::predict(downwind_lmm_fit, newdata = downwind_positive_data, random.only = TRUE)
    }

    #compute log_hatw = (1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t, for PDR (i,t) or PDR Target (i,t)
    #Depends on if there is offset term on LHS of formula(downwind_lmm_fit)
    if(length(formula.tools::lhs.vars(stats::formula(downwind_lmm_fit))) == 1){
      log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row])
    }

    if(length(formula.tools::lhs.vars(stats::formula(downwind_lmm_fit))) > 1){
      all_offset_terms = formula.tools::lhs.vars(stats::formula(downwind_lmm_fit))[2:length(formula.tools::lhs.vars(stats::formula(downwind_lmm_fit)))]
      if(length(all_offset_terms) > 1){
        log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row]) + apply(downwind_positive_data[,all_offset_terms],1,sum)
      }

      if(length(all_offset_terms) == 1){
        log_hatw = as.vector(x_z_mat[,c('(Intercept)',x_downwind_name)] %*% hatalpha_downwind + hatu[downwind_positive_useful_row]) + as.vector(downwind_positive_data[,all_offset_terms])
      }
    }


    #compute log_hatd = z_it %*% hatbeta, for PDR (i,t) or PDR Target (i,t)
    log_hatd = as.vector(x_z_mat[,setdiff(colnames(x_z_mat),c('(Intercept)',x_downwind_name))] %*% hatbeta_downwind )

    #compute log_haty_naive = (1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t + z_it %*% hatbeta, for PDR (i,t) or PDR Target (i,t)
    log_haty_naive = log_hatw + log_hatd

    #compute mu = (1/n) * sum of  (y_{i't'} / haty_{i't'}^naive ), where the sum is across PDR (i,t) or PDR Target (i,t), and n is the number of PDR (i,t) or PDR Target (i,t)  observation.
    mu = mean( y_vec / exp(log_haty_naive) )

    #compute m = Var{(1, elevation, natural_pred) %*% hatalpha_downwind + hatu_t} / Var{z_it %*% hatbeta}, where the Var is across PDR (i,t) or PDR Target (i,t).
    m = (stats::var(log_hatw))/stats::var(log_hatd)

    #compute lambda based on the formula in Ray's ISR paper
    lambda = 1 + (
      (sqrt( (1 + m)^2 + 4 * (mu - 1) * m  ) - (1 + m)) / (2 * m)
    )

    #compute hatE_it + 1 = min{2, lambda * exp(z_it %*% hatbeta)}, for PDR (i,t) or PDR Target (i,t)
    hatE_plus_one = lambda*exp(log_hatd)

    #hatR_it = y_it / (hatE_it + 1), for PDR (i,t) or PDR Target (i,t)
    hatR =  y_vec / hatE_plus_one

    #compute hatA_it = y_it - hatR_it, for PDR (i,t) or PDR Target (i,t)
    hatA = y_vec - hatR
  }

  if(attr_type == 'ThoEtAl'){
    #compute log_hatd = z_it %*% hatbeta, for PDR (i,t) or PDR Target (i,t)
    log_hatd = as.vector(x_z_mat[,setdiff(colnames(x_z_mat),c('(Intercept)',x_downwind_name))] %*% hatbeta_downwind )

    #Sigma_full is the full covariance matrix of (hatintercept, hatalpha, hatbeta)
    Sigma_full = as.matrix(stats::vcov(downwind_lmm_fit))

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








sate_est = function(downwind_positive_data, downwind_positive_target, downwind_positive_control, downwind_propensity_formula, downwind_separate_formula,
                    x_downwind_name, downwind_lmm_fit, hatalphabeta = NULL, hatu = NULL){
  #Compute different types of SATE estimate in Chambers et al. (2022)
  #When we consider hatbeta from MQ, can add another option to use either hatbeta_{0.5} or hatbeta_{hat{q}_i}
  #Note that when the instr_pred or natural_pred is used as an offset term, the sate.ipw is computed using LogRain - natural_pred, instead of LogRain only

  downwind_propensity_fit = glm(downwind_propensity_formula, family = binomial, data = downwind_positive_data)
  hatpi = predict(downwind_propensity_fit, type = "response")

  hatw_1 = (1/hatpi)/( sum( (1/hatpi) * as.numeric(downwind_propensity_fit$y) )   )
  hatw_0 = (1/(1-hatpi))/( sum( (1/ (1-hatpi)  ) * (1-as.numeric(downwind_propensity_fit$y))  )   )


  x_z_mat = stats::model.matrix(downwind_lmm_fit, data = downwind_positive_data, type = 'fixed')
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
  downwind_positive_target_lmm_fit = lme4::lmer(downwind_separate_formula, data = downwind_positive_data[downwind_positive_target,])
  downwind_positive_control_lmm_fit = lme4::lmer(downwind_separate_formula, data = downwind_positive_data[downwind_positive_control,])
  hatm_1 = predict(downwind_positive_target_lmm_fit, newdata = downwind_positive_data, re.form = NA)
  hatm_0 = predict(downwind_positive_control_lmm_fit, newdata = downwind_positive_data, re.form = NA)
  sate.ipw.ma = mean(hatm_1) - mean(hatm_0) + sum(hatw_1 * as.numeric(downwind_propensity_fit$y) * (lme4::getME(downwind_lmm_fit, 'y') - hatm_1)  ) - sum( hatw_0 * (1 - as.numeric(downwind_propensity_fit$y) ) * (lme4::getME(downwind_lmm_fit, 'y'  ) - hatm_0)  )

  #TODO: Need to check when this function is used to compute estimated sate.aipw within each bootstrap run and we are using natural_pred as offset term, should we still follow the equation (7) in JRSSA paper,
  #i.e., we replace y_i with (LogRain_i - natural_pred_i) which is captured by lme4::getME(downwind_lmm_fit, 'y') ) below that returns the response of b_downwind_lmm_fit i.e., Lograin_i^* - natural_pred_i
  sate.aipw = sum( hatw_1 * ( ( as.numeric(downwind_propensity_fit$y) * lme4::getME(downwind_lmm_fit, 'y') ) - ( (as.numeric(downwind_propensity_fit$y) - hatpi ) * hatm_1   )  )  ) - sum( hatw_0 * ( ( (1 - as.numeric(downwind_propensity_fit$y)) *  lme4::getME(downwind_lmm_fit, 'y'  )  ) -  ( (as.numeric(downwind_propensity_fit$y) - hatpi ) * hatm_0   )   )   )


  return(list(
    estimates = list(sate.mb = sate.mb,
                     sate.ipw = sate.ipw,
                     sate.ipw.l = sate.ipw.l,
                     sate.ipw.ma = sate.ipw.ma,
                     sate.aipw = sate.aipw),
    fitted_models = list(
      downwind_propensity_fit = downwind_propensity_fit,
      downwind_positive_target_lmm_fit = downwind_positive_target_lmm_fit,
      downwind_positive_control_lmm_fit = downwind_positive_control_lmm_fit
    )
  ))
}


fit_upwind_downwind_models = function(data, upwind_lmm_formula, instr_pred_name, instr_pred_type, downwind_lmm_formula, downwind_logistic_formula = NULL, upwind, downwind, positive){
  #data is the full data

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



remove_fixed_terms <- function(input_formula, vars_to_remove){
  #Function to remove some variables from a given formula - this is used to obtain the formula for fitting separate LMM to downwind_positive_target and downwind_positive_control

  # Note that this function requires downwind_lmm_formula to not use the interaction syntax such as x1*x2, but instead it should always use x1 + x2 + x1:x2


  # Flatten input_formula to a single string
  rhs_str = paste(deparse(formula.tools::rhs(input_formula)), collapse = "")
  rhs_terms = trimws(unlist(strsplit(rhs_str,'\\+')))
  rhs_terms = rhs_terms[!rhs_terms %in% vars_to_remove]

  new_rhs_str = paste(rhs_terms, collapse = ' + ')
  lhs_str = paste(deparse(formula.tools::lhs(input_formula)), collapse = "")

  # Build new input_formula
  return(stats::as.formula(paste(lhs_str, "~", new_rhs_str)))
}




