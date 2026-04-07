#' @title
#' Gauge-Day Level Observation Downwind Indicator
#'
#' @description
#' This dataset contains binary indicators of whether each gauge-day observation is downwind of each ionizer during the 2013--2018 Oman rainfall enhancement trial, as used in the analyses of Chambers et al. (2022a,b).
#'
#'
#' @format A binary matrix with 122259 rows (corresponding to gauge-days in \code{\link{oman}} dataset) and 10 variables:
#' \describe{
#'    \item{H1 -- H10}{Binary indicators for each ionizer, showing whether a gauge-day observation is downwind of that ionizer: 1 = downwind; 0 = not downwind; NA = ionizer not yet deployed}
#'}
#'
#'@references
#'\itemize{
#'  \item Chambers, R., Beare, S., Peak, S. and Al-Kalbani, M. (2022a) Nudging a Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall Enhancement Trial in Oman. \emph{International Statistical Review}, 90: 346–373.
#'  \item Chambers, R., Ranjbar, S., Salvati, N., and Pacini, B. (2022b) Weighting, Informativeness and Causal Inference, with an Application to Rainfall Enhancement. \emph{Journal of the Royal Statistical Society Series A: Statistics in Society}, 185: 1584–1612
#'}
#'
"gaugeday_downwind"
