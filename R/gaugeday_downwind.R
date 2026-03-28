#' @title
#' Gauge-Day Level Observation Downwind Indicator
#'
#' @description
#' This dataset contains binary indicators of whether each gauge-day observation is downwind of each ionizer during the 2013--2018 Oman rainfall enhancement trial, as used in the analyses of Chambers et al. (2022a,b).
#'
#'
#' @format A binary matrix with 122259 rows (gauge-days) and 10 variables:
#' \describe{
#'    \item{H1 -- H10}{Binary indicators for each ionizer, showing whether a gauge-day observation is downwind of that ionizer: 1 = downwind; 0 = not downwind; NA = ionizer not yet deployed}
#'
#'
#'
#'
"gaugeday_downwind"
