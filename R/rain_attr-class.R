#TODO: Modify the following documentation to match what is actually being done

#' @title Class "rain_attr" of Rainfall Attribution Analysis
#' @description
#' Objects of class `rain_attr` represent results from the two-stage LMM-based
#' rainfall enhancement analysis, including fitted models, bootstrap/permutation
#' summaries, and diagnostic tools.
#'
#' @details
#' A `rain_attr` object contains fitted mixed models (LMMs), GLMs, and associated
#' results such as attribution estimates, SATE, and variance components.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{summary(object, ...)}}{Summarizes attribution, SATE, and model details.}
#'   \item{\code{print(object, ...)}}{Prints a short description of the object.}
#'   \item{\code{plot(object, ...)}}{Diagnostic and bootstrap/permutation plots.}
#'   \item{\code{residuals(object, ...)}}{Extracts residuals for a given model.}
#'   \item{\code{predict(object, ...)}}{Predicts fitted values from the models.}
#'   \item{\code{coef(object, ...)}}{Extracts fixed effects or coefficients.}
#'   \item{\code{varcomp(object, ...)}}{Extracts variance components for LMMs.}
#' }
#'
#' @seealso \code{\link{rain_attr}}, \code{\link{summary.rain_attr}}
#'
#' @docType class
#' @name rain_attr-class
NULL
