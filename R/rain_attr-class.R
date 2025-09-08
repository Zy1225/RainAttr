#TODO: Modify the following documentation to match what is actually being done
#TODO: Continue adding description of each @param

#' @title Class "rain_attr" of Two-Stage LMM Fitted to Rainfall Enhancement Trial Data
#' @description
#' Objects of class `rain_attr` represent results from the two-stage LMM-based
#' rainfall enhancement analysis, created by calls to \code{\link{rain_attr}}.
#'
#' @details
#' A `rain_attr` object contains fitted LMMs:
#' \describe{
#'  \item{\code{upwind_lmm}}{Upwind (first stage) LMM.}
#'  \item{\code{downwind_lmm}}{Downwind (second stage) LMM.}
#'  \item{\code{downwind_target_lmm}}{Downwind (second stage) treatment-only LMM.}
#'  \item{\code{downwind_control_lmm}}{Downwind (second stage) control-only LMM.}
#' }
#' and fitted GLMs:
#' \describe{
#'  \item{\code{downwind_logistic}}{Downwind (second stage) logistic model of rainfall event indicator.}
#'  \item{\code{downwind_propensity}}{Downwind (second stage) propensity score model for the treatment indicator.}
#' }
#'
#' The object also includes attribution estimates, SATE, and results from bootstrap and permutation-based procedures if \code{bootstrap = TRUE}
#' or \code{permutation = TRUE} were specified in the original call.
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
#'
#' @param object An \bold{R} object of class \code{\link{rain_attr-class}}, i.e., as resulting from \code{\link{rain_attr}()}.
#' @param model A character string specifying which fitted model to focus on. Must be one of \code{upwind_lmm}, \code{downwind_lmm}, \code{downwind_target_lmm}, \code{downwind_control_lmm}, \code{downwind_logistic}, or \code{downwind_propensity}. See "Details" for more information.
#'
#' @seealso \code{\link{rain_attr}}, \code{\link{summary.rain_attr}}
#'
#' @docType class
#' @name rain_attr-class
NULL
