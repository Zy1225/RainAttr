#' @title
#' Daily Oman Ionizer Operation Schedule
#'
#' @description
#' This dataset contains the daily operation schedule of the ionizers during the 2013--2018 Oman rainfall enhancement trial, as used in the analyses of Chambers et al. (2022a,b).
#'
#'
#' @format A data frame with 740 rows (days) and 12 variables:
#' \describe{
#'    \item{TrialDay}{Identifier for each day in the six-year trial}
#'    \item{Year}{Character variable indicating the year}
#'    \item{H1 -- H10}{Binary indicator variables for the operating status of the ionizers: 1 = ionizer active on that day, 0 = ionizer inactive, NA = ionizer not yet deployed}
#'
#' }
#'
#'
#'
#'
"ionizer_operation"
