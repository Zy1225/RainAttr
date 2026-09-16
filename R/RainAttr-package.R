#' RainAttr: Attribution and Sample Average Treatment Effect for Rainfall Enhancement Trial Data
#'
#' The RainAttr package provides tools for estimation and inference of attribution and sample
#' average treatment effect in rainfall enhancement trials using a two-stage linear mixed model approach
#' based on Chambers et al. (2022).
#'
#' The main functions are:
#'
#' \itemize{
#'   \item \code{\link{rain_attr}} for fitting the two-stage model,
#'   estimating attribution and SATEs, and optionally performing bootstrap-based and permutation-based inference.
#'   \item \code{\link{bootstrap_opt}} for configuring bootstrap inference.
#'   \item \code{\link{permutation_opt}} for configuring permutation inference.
#'   \item \code{\link{eda}} for exploratory data analysis.
#' }
#'
#' See the package vignette for a complete worked example based on the Oman
#' rainfall enhancement trial.
#'
#'@references
#'\itemize{
#'  \item Chambers, R., Beare, S., Peak, S. and Al-Kalbani, M. (2022) Nudging a Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall Enhancement Trial in Oman. \emph{International Statistical Review}, 90: 346–373.
#'}
"_PACKAGE"
