
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
#' When \code{discretize_rain = TRUE}, bootstrap samples that satisfy \eqn{ Rain_{ij}^* \in (0, 0.3] } are replaced by 0.2, \eqn{ Rain_{ij}^* \in (0.3, 0.5] } are replaced by 0.4, \eqn{ Rain_{ij}^* \in (0.5, 0.7] } are replaced by 0.6, and \eqn{ Rain_{ij}^* \in (0.7, 0.9] } are replaced by 0.8.
#'
#' When \code{winsorize_individual_rain = TRUE}, bootstrap samples that satisfy \eqn{Rain_{ij}^* >  } \code{individual_rain_interval[2]} are replaced by random numbers drawn from a uniform distribution over the interval \eqn{[}\code{individual_rain_interval[1]}, \code{individual_rain_interval[2]}\eqn{]}.
#'
#' When \code{winsorize_total_rain = TRUE}, if \eqn{ \sum_{(i,j)} Rain_{ij}^* \notin [} \code{total_rain_interval[1]}, \code{total_rain_interval[2]} \eqn{]} where the summation is over the subset of observations in \code{ori_data} satisfying \code{downwind} and \eqn{L_{ij}^* = 1}, then all bootstrap samples of \eqn{Rain_{ij}^*} are rescaled by a common factor of \code{runif(n = 1, min = total_rain_interval[1], max = total_rain_interval[2])} \eqn{ / \sum_{(i,j)} Rain_{ij}^*  }.
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
#'   (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param bootstrap_type A character string specifying the type of bootstrap.
#'   Must be one of \code{"REB0"}, \code{"REB1"}, \code{"REB2"}, \code{"PREB0"},
#'   \code{"PREB1"}, \code{"PREB2"}, or \code{"MREB1"}.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param bootstrap_zero Logical. If \code{TRUE}, the optional first-level bootstrap is performed to generate bootstrap samples of binary rainfall event indicators.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param positive_prob_threshold An optional numeric value between 0 and 1 specifying the probability threshold for generating bootstrap samples of binary rainfall event indicators. Probabilities below this threshold are set to zero.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param discretize_rain Logical. If \code{TRUE}, rainfall values are discretized in bootstrap resamples.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param winsorize_individual_rain Logical. If \code{TRUE}, individual rainfall values in bootstrap samples that exceed the upper bound specified by \code{individual_rain_interval} are replaced with random draws from a uniform distribution over \code{[individual_rain_interval[1], individual_rain_interval[2]]}.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param individual_rain_interval Numeric vector of length 2 specifying the lower and upper bounds for adjusting bootstrapped individual rainfall values that are too large when \code{winsorize_individual_rain = TRUE}.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param winsorize_total_rain Logical. If \code{TRUE}, all individual rainfall values in each bootstrap sample are proportionally rescaled so that the total equals a random number drawn uniformly from
#'   \code{[total_rain_interval[1], total_rain_interval[2]]} whenever the total bootstrapped rainfall falls outside this interval.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param total_rain_interval Numeric vector of length 2 specifying the lower and upper bounds for adjusting the total of bootstrapped rainfall values when \code{winsorize_total_rain = TRUE}.
#'   (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param bootstrap_seed An integer specifying the random seed for the bootstrap procedure. Reproducibility is guaranteed only if \code{bootstrap_parallel} is the same, since parallel execution changes the order of random number generation.
#' (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param bootstrap_parallel Logical. If \code{TRUE}, each bootstrap run is executed in parallel across multiple workers. If \code{FALSE}, they are run sequentially.
#' (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
#' @param bootstrap_parallel_num_worker An integer specifying the number of parallel workers to use when \code{bootstrap_parallel = TRUE}.
#' (User-configurable bootstrap option using \code{\link{bootstrap_opt}})
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
#' @param attr_type A character string specifying the type of attribution estimates. Must be one of \code{"ChambersEtAl"}, \code{"ChambersEtAl_No_Winsorize"}, \code{"ThoEtAl"}, or \code{"No"}. See \code{\link{rain_attr}} for more information.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param x_downwind_name A character vector containing variable names from the right hand side of \code{downwind_lmm_formula}, for those variables that are not related to ionizers (treatment). The intercept is always included and does not need to be specified.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param target_only Logical. If \code{TRUE} the attribution estimates are computed based on only target observations. If \code{FALSE} the attribution estimates are computed based on both treatment and control observations.
#'   (Internal argument set automatically when using \code{\link{rain_attr}})
#' @param downwind_propensity_formula A two sided linear formula object to be used in \code{\link{glm}} with \code{family = "binomial"}, for fitting a propensity score model to the treatment indicators of downwind (second stage) observations.
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
#' @seealso \code{\link{rain_attr}} for the main function, \code{\link{bootstrap_opt}} for specifying bootstrap options
#'
#' @references
#'\itemize{
#'  \item Chambers, R. and Chandra, H. (2013). A Random Effect Block Bootstrap for Clustered Data. \emph{Journal of Computational and Graphical Statistics}, 22, 452–470.
#'  \item Tho, Z. Y., Chambers, R., and Welsh, A. H. (2025) Adjusted Random Effect Block Bootstraps for Highly Unbalanced Clustered Data.
#'     \href{https://arxiv.org/abs/2510.07770}{arXiv:2510.07770}.
#'}
#'


#TODO: Consider to add another variation of 'PREB2' and 'REB2' for adjusting downwind_lmm_fit's fixef as well as random effects, and plug these corrected estimates to compute hatattr and hatsate, rather than directly centering hatattr and hatsate
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


  r_vec = lme4::getME(ori_fitted_models$downwind_lmm_fit, 'y')  - predict(ori_fitted_models$downwind_lmm_fit, re.form = NA)
  group_name = names(lme4::getME(ori_fitted_models$downwind_lmm_fit, "flist"))
  ori_downwind_positive_group = ori_data[downwind & ori_positive , group_name]
  ori_downwind_positive_group_label = unique(ori_downwind_positive_group)
  ori_D_groups =  length(ori_downwind_positive_group_label)

  if(bootstrap_type %in% c('REB0', 'REB2', 'PREB0', 'PREB2')){
    hat.u = sapply(ori_downwind_positive_group_label, function(x){
      mean(r_vec[ori_downwind_positive_group == x])
    })


    hat.e = rep(0,sum(downwind & ori_positive))
    for(h in 1:ori_D_groups){
      hat.e[ori_downwind_positive_group==ori_downwind_positive_group_label[h]] = r_vec[ori_downwind_positive_group==ori_downwind_positive_group_label[h]]- hat.u[h]
    }

    final.hat.u = hat.u
    final.hat.e = hat.e
  }

  if(bootstrap_type %in%  c('REB1','PREB1', 'MREB1')){
    hat.u = sapply(ori_downwind_positive_group_label, function(x){
      mean(r_vec[ori_downwind_positive_group == x])
    })

    hat.e = rep(0,sum(downwind & ori_positive))
    for(h in 1:ori_D_groups){
      hat.e[ori_downwind_positive_group==ori_downwind_positive_group_label[h]] = r_vec[ori_downwind_positive_group==ori_downwind_positive_group_label[h]]- hat.u[h]
    }


    vc_df = as.data.frame(lme4::VarCorr(ori_fitted_models$downwind_lmm_fit))
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


  if(!bootstrap_parallel){
    if(!is.null(bootstrap_seed)){
      set.seed(bootstrap_seed)
    }
    for(b in 1:B_bootstrap){
      tryCatch({

        if(bootstrap_zero){
          b_downwind_positive = (runif(num_downwind, min = 0, max = 1) < downwind_positive_prob)
          b_downwind_logistic_fit = glm(b_downwind_positive ~ model.matrix(ori_fitted_models$downwind_logistic_fit$formula, data = ori_data[downwind,]) - 1, family = 'binomial')
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

        b_y = rep(NA, b_num_downwind_positive)

        b_donor_group_label = sample(x = ori_downwind_positive_group_label,
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


          b_y[target.units] = b_fitted[target.units] + b_u[h] + final.hat.e[donating.units]
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
          b_raw_y[b_raw_y<0.3] = 0.2
          b_raw_y[(b_raw_y>0.3)&(b_raw_y<0.5)] = 0.4
          b_raw_y[(b_raw_y>0.5)&(b_raw_y<0.7)] = 0.6
          b_raw_y[(b_raw_y>0.7)&(b_raw_y<0.9)] = 0.8
        }



        if(winsorize_individual_rain){
          b_raw_y[b_raw_y> individual_rain_interval[2]] = individual_rain_interval[1] + (individual_rain_interval[2] - individual_rain_interval[1]) * runif(n=sum(b_raw_y> individual_rain_interval[2]))
        }

        if(winsorize_total_rain){
          if(sum(b_raw_y)< total_rain_interval[1] | sum(b_raw_y)> total_rain_interval[2]){
            b_raw_y = b_raw_y*(runif(n=1,min=total_rain_interval[1], max=total_rain_interval[2]))/sum(b_raw_y)
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
    parallel::clusterExport(cl, varlist = c('attr_est', 'sate_est'), envir = asNamespace("RainAttr"))
    doParallel::registerDoParallel(cl)
    if(!is.null(bootstrap_seed)){
      doRNG::registerDoRNG(seed = bootstrap_seed)
    }
    `%dopar%` = foreach::`%dopar%`

    results = foreach::foreach(b = 1:B_bootstrap, .packages = c("lme4","rlang","formula.tools"), .errorhandling = 'remove') %dopar% {
      if(bootstrap_zero){
        b_downwind_positive = (runif(num_downwind, min = 0, max = 1) < downwind_positive_prob)
        b_downwind_logistic_fit = glm(b_downwind_positive ~ model.matrix(ori_fitted_models$downwind_logistic_fit$formula, data = ori_data[downwind,]) - 1, family = 'binomial')
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

      b_y = rep(NA, b_num_downwind_positive)

      b_donor_group_label = sample(x = ori_downwind_positive_group_label,
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


        b_y[target.units] = b_fitted[target.units] + b_u[h] + final.hat.e[donating.units]
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
        b_raw_y[b_raw_y<0.3] = 0.2
        b_raw_y[(b_raw_y>0.3)&(b_raw_y<0.5)] = 0.4
        b_raw_y[(b_raw_y>0.5)&(b_raw_y<0.7)] = 0.6
        b_raw_y[(b_raw_y>0.7)&(b_raw_y<0.9)] = 0.8
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
#' @param bootstrap_seed An integer specifying the random seed for the bootstrap procedure. Reproducibility is guaranteed only if \code{bootstrap_parallel} is the same, since parallel execution changes the order of random number generation. Default is \code{NULL}, meaning no seed is set internally and users should call \code{set.seed()} beforehand to ensure reproducibility.
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
#' # Specifically: bootstrap_type = 'PREB1' as proposed by Tho et al. (2025) Adjusted Random Effect Block Bootstraps for Highly Unbalanced Clustered Data. arXiv:2510.07770.
#' boot_options = bootstrap_opt()
#' str(boot_options)
#'
#' #Bootstrap option with parallelization over (parallel::detectCores() - 1) number of workers and seed = 1 for reproducibility
#' boot_options_parallel = bootstrap_opt(
#'   bootstrap_seed = 1,
#'   bootstrap_parallel = TRUE,
#'   bootstrap_parallel_num_worker = parallel::detectCores() - 1
#' )
#' str(boot_options_parallel)
#' @export

bootstrap_opt = function(B_bootstrap = 10000,
                         bootstrap_type = 'PREB1',
                         bootstrap_zero = T,
                         positive_prob_threshold = NULL,
                         discretize_rain = T,
                         winsorize_individual_rain = T,
                         individual_rain_interval = c(100,175),
                         winsorize_total_rain = T,
                         total_rain_interval = c(6000,60000),
                         bootstrap_seed = NULL,
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
