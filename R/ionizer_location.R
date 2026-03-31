#' @title
#' Oman Ionizer Location and Deployment Timing
#'
#' @description
#' This dataset contains the latitude and longitude of the ten ionizers and their deployment timing during the 2013--2018 Oman rainfall enhancement trial, as used in the analyses of Chambers et al. (2022a,b).
#'
#'
#' @format A data frame with 10 rows (ionizers) and 9 variables:
#' \describe{
#'    \item{Ionizer}{Identifier for each ionizer}
#'    \item{Latitude, Longitude}{Geographic coordinates of each ionizer}
#'    \item{2013 -- 2018}{Binary indicator variables for the deployment status of each ionizer in each year: 1 = deployed, 0 = not yet deployed}
#'
#' }
#'
#'@references
#'\itemize{
#'  \item Chambers, R., Beare, S., Peak, S. and Al-Kalbani, M. (2022a) Nudging a Pseudo-Science Towards a Science—The Role of Statistics in a Rainfall Enhancement Trial in Oman. \emph{International Statistical Review}, 90: 346–373.
#'  \item Chambers, R., Ranjbar, S., Salvati, N., and Pacini, B. (2022b) Weighting, Informativeness and Causal Inference, with an Application to Rainfall Enhancement. \emph{Journal of the Royal Statistical Society Series A: Statistics in Society}, 185: 1584–1612
#'}
#'
#'
"ionizer_location"
