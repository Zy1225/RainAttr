#' @rdname rain_attr-class
#' @export
coef.rain_attr = function(object, model = "downwind_lmm", ...){
  # Map model names to stored fitted models
  model_map = list(
    upwind_lmm = object$all_fitted_models$upwind_lmm_fit,
    downwind_lmm = object$all_fitted_models$downwind_lmm_fit,
    downwind_logistic = object$all_fitted_models$downwind_logistic_fit,
    downwind_propensity = object$all_fitted_models$downwind_propensity_fit,
    downwind_target_lmm = object$all_fitted_models$downwind_positive_target_lmm_fit,
    downwind_control_lmm = object$all_fitted_models$downwind_positive_control_lmm_fit
  )

  if (!model %in% names(model_map)) {
    stop(paste("Invalid model. Choose one of:", paste(names(model_map), collapse = ", ")))
  }

  if(model == 'downwind_logistic'){
    if (is.null(object$all_fitted_models$downwind_logistic_fit)) {
      stop("No downwind_logistic_formula has been specified in the original call of rain_attr.
               Please rerun rain_attr with downwind_logistic_formula specified.")
    }
  }

  selected_model = model_map[[model]]

  if(model %in%  c('upwind_lmm', 'downwind_lmm', 'downwind_target_lmm', 'downwind_control_lmm' )){
    output = list(
      fixef_coef = lme4::fixef(selected_model),
      ranef_coef = lme4::ranef(selected_model)[[1]]
    )
  }else{
    output = coef(selected_model)
  }

  return(output)
}


#' @rdname rain_attr-class
#' @export
print.rain_attr <- function(object, ...) {
  cat("Two Stage LMM Rainfall Enhancement Analysis Result\n")
  cat("======================================================================\n\n")

  # Point estimates
  cat("Point Estimates:\n")

  # Attribution (hatattr) as %
  hatattr_values <- c(
    apo = if(!is.null(object$hatattr$apo)) paste0(round(object$hatattr$apo * 100, 2), "%") else NA,
    apl = if(!is.null(object$hatattr$apl)) paste0(round(object$hatattr$apl * 100, 2), "%") else NA
  )
  cat("Attribution (%) Assuming Log-Rainfall being Modelled:\n")
  print(noquote(hatattr_values))
  cat("\n")

  # SATE estimates (hatsate)
  sate_names <- c("sate.mb", "sate.ipw", "sate.ipw.l", "sate.ipw.ma", "sate.aipw")
  sate_values <- sapply(sate_names, function(nm) if(!is.null(object$hatsate[[nm]])) object$hatsate[[nm]] else NA)

  cat("SATE Estimates (hatsate):\n")
  print(noquote(format(sate_values, digits = 4)))

  cat("======================================================================\n\n")
  # Indicate whether bootstrap or permutation was carried out
  cat("Inference:\n")
  if (!is.null(object$bootstrap_result)) {
    cat("Bootstrap inference has been carried out. See summary() for detailed results.\n\n")
  } else {
    cat("Bootstrap inference has NOT been carried out.\n\n")
  }

  if (!is.null(object$permutation_result)) {
    cat("Permutation inference has been carried out. See summary() for detailed results.\n\n")
  } else {
    cat("Permutation inference has NOT been carried out.\n\n")
  }

  cat("======================================================================\n\n")
  # Upwind formula and stats
  upwind_formula <- object$args$upwind_lmm_formula
  upwind_fit <- object$all_fitted_models$upwind_lmm_fit

  cat("Upwind (First Stage) LMM Formula:\n")
  print(upwind_formula)

  cat('\n')
  cat("Data subset used: ")
  cat(deparse(object$args$data),"[", deparse(object$args$upwind_subset), " & ", deparse(object$args$positive_subset), ", ]\n")


  n_obs_up <- nobs(upwind_fit)
  n_groups_up <- length(unique(lme4::getME(upwind_fit, "flist")[[1]]))
  cat(sprintf("Number of observations: %d, ", n_obs_up))
  cat(sprintf("Number of groups: %d\n\n", n_groups_up))


  # Random effects
  cat("Random effects:\n")
  re_up <- lme4::VarCorr(upwind_fit)
  print(re_up)
  cat("\n")

  # Fixed effects
  cat("Fixed effects:\n")
  fe_up <- lme4::fixef(upwind_fit)
  print(round(fe_up, 4))
  cat("\n======================================================================\n\n")

  # Downwind formula and stats
  downwind_formula <- object$args$downwind_lmm_formula
  downwind_fit <- object$all_fitted_models$downwind_lmm_fit

  cat("Downwind (Second Stage) LMM Formula:\n")
  print(downwind_formula)

  cat('\n')
  cat("Data subset used: ")
  cat(deparse(object$args$data),"[", deparse(object$args$downwind_subset), " & ", deparse(object$args$positive_subset), ", ]\n")


  n_obs_down <- nobs(downwind_fit)
  n_groups_down <- length(unique(lme4::getME(downwind_fit, "flist")[[1]]))
  cat(sprintf("Number of observations: %d, ", n_obs_down))
  cat(sprintf("Number of groups: %d\n\n", n_groups_down))

  # Random effects
  cat("Random effects:\n")
  re_down <- lme4::VarCorr(downwind_fit)
  print(re_down)
  cat("\n")

  # Fixed effects
  cat("Fixed effects:\n")
  fe_down <- lme4::fixef(downwind_fit)
  print(round(fe_down, 4))

  invisible(object)
}



#' @rdname rain_attr-class
#' @export

residuals.rain_attr <- function(object, model = "downwind_lmm",
                                residual_type = NULL, residual_scaled = TRUE, ...) {
  # Match arguments
  model_map <- list(
    upwind_lmm = object$all_fitted_models$upwind_lmm_fit,
    downwind_lmm = object$all_fitted_models$downwind_lmm_fit,
    downwind_logistic = object$all_fitted_models$downwind_logistic_fit,
    downwind_propensity = object$all_fitted_models$downwind_propensity_fit,
    downwind_target_lmm = object$all_fitted_models$downwind_positive_target_lmm_fit,
    downwind_control_lmm = object$all_fitted_models$downwind_positive_control_lmm_fit
  )

  if (!model %in% names(model_map)) {
    stop(paste("Invalid model. Choose one of:", paste(names(model_map), collapse = ", ")))
  }

  if(model == 'downwind_logistic'){
    if (is.null(object$all_fitted_models$downwind_logistic_fit)) {
      stop("No downwind_logistic_formula has been specified in the original call of rain_attr.
               Please rerun rain_attr with downwind_logistic_formula specified.")
    }
  }

  selected_model <- model_map[[model]]


  if (inherits(selected_model, "merMod")) {
    if (is.null(residual_type)) residual_type <- "response"  # Default for LMM
    valid_types <- c("working", "response", "deviance", "pearson")
    if (!residual_type %in% valid_types) {
      stop(paste("For LMM models, residual_type must be one of:", paste(valid_types, collapse = ", ")))
    }
    res <- resid(selected_model, type = residual_type, scaled = residual_scaled)
  }

  if (inherits(selected_model, "glm")) {
    # GLM models
    if (is.null(residual_type)) residual_type <- "deviance"  # Default for GLM
    valid_types <- c("deviance", "pearson", "response", "working", "partial")
    if (!residual_type %in% valid_types) {
      stop(paste("For GLM models, residual_type must be one of:", paste(valid_types, collapse = ", ")))
    }
    res <- residuals(selected_model, type = residual_type)

  }

  return(res)
}


#' @rdname rain_attr-class
#' @export

fitted.rain_attr = function(object, model = "downwind_lmm", ...){
  return(
    predict.rain_attr(object, model = model, predict_type = "response")
  )
}




#' Generic for variance components
#' @rdname rain_attr-class
#' @usage NULL
#' @export
varcomp <- function(object, ...) {
  UseMethod("varcomp")
}


#' @rdname rain_attr-class
#' @export
varcomp.rain_attr <- function(object, ...) {

  # Extract LMMs
  lmm_list <- list(
    upwind_lmm = object$all_fitted_models$upwind_lmm_fit,
    downwind_lmm = object$all_fitted_models$downwind_lmm_fit,
    downwind_target_lmm = object$all_fitted_models$downwind_positive_target_lmm_fit,
    downwind_control_lmm = object$all_fitted_models$downwind_positive_control_lmm_fit
  )


  res_mat <- matrix(NA, nrow = 4, ncol = 2)
  rownames(res_mat) <- names(lmm_list)


  # Fill in variance components
  for (i in seq_along(lmm_list)) {
    res_mat[i,] = as.data.frame(lme4::VarCorr(lmm_list[[i]]))[,'vcov']
  }
  colnames(res_mat) <- as.data.frame(lme4::VarCorr(lmm_list[[i]]))[,'grp']

  return(res_mat)
}



#' @rdname rain_attr-class
#' @export
#'
predict.rain_attr = function(object, newdata = NULL, model = "downwind_lmm",
                             re_include = TRUE, fixef_include = TRUE, allow.new.levels = FALSE,
                             predict_type = "link", ...) {
  # Map model names to stored fitted models
  model_map = list(
    upwind_lmm = object$all_fitted_models$upwind_lmm_fit,
    downwind_lmm = object$all_fitted_models$downwind_lmm_fit,
    downwind_logistic = object$all_fitted_models$downwind_logistic_fit,
    downwind_propensity = object$all_fitted_models$downwind_propensity_fit,
    downwind_target_lmm = object$all_fitted_models$downwind_positive_target_lmm_fit,
    downwind_control_lmm = object$all_fitted_models$downwind_positive_control_lmm_fit
  )

  if (!model %in% names(model_map)) {
    stop(paste("Invalid model. Choose one of:", paste(names(model_map), collapse = ", ")))
  }

  if(model == 'downwind_logistic'){
    if (is.null(object$all_fitted_models$downwind_logistic_fit)) {
      stop("No downwind_logistic_formula has been specified in the original call of rain_attr.
               Please rerun rain_attr with downwind_logistic_formula specified.")
    }
  }

  selected_model = model_map[[model]]

  # Predictions for lme4::merMod
  if (inherits(selected_model, "merMod")) {
    if(re_include){
      fit = predict(selected_model, newdata = newdata,
                    re.form = NULL,
                    random.only = !fixef_include,
                    allow.new.levels = allow.new.levels, ...)
    }else{
      fit = predict(selected_model, newdata = newdata,
                    re.form = NA,
                    random.only = !fixef_include,
                    allow.new.levels = allow.new.levels, ...)
    }
  }

  # Predictions for glm
  if (inherits(selected_model, "glm")) {
    if(!predict_type %in% c("link", "response", "terms")){
      stop(paste("Invalid predict_type. Choose one of:", paste(c("link", "response", "terms"), collapse = ", ")))
    }

    fit = predict(selected_model, newdata = newdata,
                  type = predict_type, ...)
  }

  return(fit)
}

#' @rdname rain_attr-class
#' @export

plot.rain_attr = function(object, plot_type = c("bootstrap", "permutation"), plot_quantity = c("attr", "sate"),
                          model = 'downwind_lmm', residual_type = NULL, residual_scaled = TRUE,
                          re_include = TRUE, fixef_include = TRUE, allow.new.levels = FALSE, predict_type = "link",
                          ...) {

  allowed_single_type <- c("model", "bootstrap", "permutation")

  if (length(plot_type) == 1) {
    if (!plot_type %in% allowed_single_type) {
      stop(sprintf("Invalid plot_type. Must be one of: %s",
                   paste(allowed_single_type, collapse = ", ")))
    }
  } else if (length(plot_type) == 2) {
    if (!setequal(plot_type, c("bootstrap", "permutation"))) {
      stop("When plot_type has length 2, it must be exactly c('bootstrap','permutation').")
    }
  } else {
    stop("plot_type must be either length 1 (model, bootstrap, or permutation) or length 2 (bootstrap, permutation).")
  }


  allowed_single_quant <- c("attr", "sate")

  if (length(plot_quantity) == 1) {
    if (!plot_quantity %in% allowed_single_quant) {
      stop(sprintf("Invalid plot_quantity. Must be one of: %s",
                   paste(allowed_single_quant, collapse = ", ")))
    }
  } else if (length(plot_quantity) == 2) {
    if (!setequal(plot_quantity, c("attr", "sate"))) {
      stop("When plot_quantity has length 2, it must be exactly c('attr','sate').")
    }
  } else {
    stop("plot_quantity must be either length 1 (attr or sate) or length 2 (attr, sate).")
  }


  if (any(plot_type %in% c("bootstrap", "permutation"))) {

    plot_lists = list()
    col_labels = c()

    for (pt in c("bootstrap", "permutation")) {

      if (!(pt %in% plot_type)) {
        next
      }

      if (pt == "bootstrap"){
        plot_obj = object$bootstrap_plot_result
      }else{
        plot_obj = object$permutation_plot_result
      }

      if (is.null(plot_obj)) {
        stop(sprintf("No %s has been carried out. Please rerun rain_attr with %s = TRUE.",
                     pt, pt))
      }

      pt_plots = list()

      for (pq in plot_quantity) {

        if (pq == "attr") {
          plots_to_show = plot_obj$hatattr
        }  else {
          plots_to_show = plot_obj$hatsate
        }

        pt_plots = c(pt_plots, plots_to_show)
      }

      plot_lists[[pt]] = pt_plots

      # Capitalize first letter for column label
      col_labels = c(col_labels, paste0(toupper(substr(pt, 1, 1)), substr(pt, 2, nchar(pt))))
    }

    ncol = length(plot_lists)
    nrow = max(sapply(plot_lists, length))

    # Arrange each column vertically
    column_arranges = lapply(names(plot_lists), function(pt) {
      ggpubr::ggarrange(
        plotlist = plot_lists[[pt]],
        ncol = 1,
        nrow = length(plot_lists[[pt]]),
        labels = NULL
      )
    })

    combined_plot = do.call(
      ggpubr::ggarrange,
      c(column_arranges,
        list(ncol = ncol, nrow = 1, labels = col_labels))
    )


    print(combined_plot)
  }

  if ('model' %in% plot_type) {

    model_map = list(
      upwind_lmm = object$all_fitted_models$upwind_lmm_fit,
      downwind_lmm = object$all_fitted_models$downwind_lmm_fit,
      downwind_logistic = object$all_fitted_models$downwind_logistic_fit,
      downwind_propensity = object$all_fitted_models$downwind_propensity_fit,
      downwind_target_lmm = object$all_fitted_models$downwind_positive_target_lmm_fit,
      downwind_control_lmm = object$all_fitted_models$downwind_positive_control_lmm_fit
    )

    plot_title_map = list(
      upwind_lmm = 'Upwind LMM',
      downwind_lmm = 'Downwind LMM',
      downwind_logistic = 'Downwind Logistic Model',
      downwind_propensity = 'Downwind Propensity Score Model',
      downwind_target_lmm = 'Downwind Target-Only LMM',
      downwind_control_lmm = 'Downwind Control-Only LMM'
    )

    if (!model %in% names(model_map)) {
      stop(paste("Invalid model. Choose one of:", paste(names(model_map), collapse = ", ")))
    }

    if(model == 'downwind_logistic'){
      if (is.null(object$all_fitted_models$downwind_logistic_fit)) {
        stop("No downwind_logistic_formula has been specified in the original call of rain_attr.
               Please rerun rain_attr with downwind_logistic_formula specified.")
      }
    }

    selected_model = model_map[[model]]
    plot_title = plot_title_map[[model]]


    # Extract residuals and fitted values
    res = residuals(object, model = model, residual_type = residual_type, residual_scaled = residual_scaled)
    fit = predict(object, newdata = NULL, model = model, re_include = re_include, fixef_include = fixef_include, allow.new.levels = allow.new.levels, predict_type = predict_type)



    #
    if (model == "downwind_lmm") {
      downwind_positive_target = rlang::eval_tidy(
        object$args$downwind_target_subset,
        data = object$data[rlang::eval_tidy(object$args$downwind_subset, data = object$data) &
                             rlang::eval_tidy(object$args$positive_subset, data = object$data), ]
      )

      control_label = as.character(object$args$downwind_control_subset)[3]
      target_label = as.character(object$args$downwind_target_subset)[3]

      df = data.frame(
        Fitted = fit,
        Residuals = res,
        Group = factor(ifelse(downwind_positive_target, target_label, control_label))
      )



      legend_title = as.character(object$args$downwind_target_subset)[2]


      color_vals <- c(
        setNames("#1f77b4", target_label),
        setNames("#ff7f0e", control_label)
      )

      shape_vals <- c(
        setNames(16, target_label),
        setNames(17, control_label)
      )

      # Residuals vs Fitted
      p1 = ggplot2::ggplot(df, ggplot2::aes(x = Fitted, y = Residuals, color = Group, shape = Group)) +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste("Residuals vs Fitted for", plot_title),
                      y = paste0("Residuals (type=", ifelse(is.null(residual_type), 'response', residual_type), ", scaled=", residual_scaled, ")"),
                      x = paste0("Fitted (fixef_include =", fixef_include, ", re_include =", re_include, ")"),
                      color = legend_title, shape = legend_title) +
        ggplot2::theme_bw() +
        ggplot2::scale_color_manual(values = color_vals) +
        ggplot2::scale_shape_manual(values = shape_vals) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 9),
                       axis.title.x = ggplot2::element_text(size = 8),
                       axis.title.y = ggplot2::element_text(size = 8))

      # QQ plot
      qq_vals = qqnorm(res, plot.it = FALSE)
      qq_df = data.frame(Theoretical = qq_vals$x, Sample = qq_vals$y, Group = df$Group)

      q_sample = quantile(res, probs = c(0.25, 0.75))
      q_theory = qnorm(c(0.25, 0.75))
      qqline_slope = diff(q_sample) / diff(q_theory)
      qqline_intercept = q_sample[1] - qqline_slope * q_theory[1]


      p2 = ggplot2::ggplot(qq_df, ggplot2::aes(x = Theoretical, y = Sample, color = Group, shape = Group)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(intercept = qqline_intercept, slope = qqline_slope, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste("Normal Q-Q Plot for", plot_title, "Residuals"),
                      color = legend_title, shape = legend_title) +
        ggplot2::theme_bw() +
        ggplot2::scale_color_manual(values = color_vals) +
        ggplot2::scale_shape_manual(values = shape_vals) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 9),
                       axis.title.x = ggplot2::element_text(size = 8),
                       axis.title.y = ggplot2::element_text(size = 8))

      # Arrange side by side with shared legend
      print(ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom"))

    }

    if(model %in% c('upwind_lmm', 'downwind_target_lmm', 'downwind_control_lmm')){
      df = data.frame(Fitted = fit, Residuals = res)


      #Residuals vs Fitted
      p1 = ggplot2::ggplot(df, ggplot2::aes(x = Fitted, y = Residuals)) +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste("Residuals vs Fitted for", plot_title),
                      y = paste0("Residuals (type=", ifelse(is.null(residual_type), 'response', residual_type), ", scaled=", residual_scaled, ")"),
                      x = paste0("Fitted (fixef_include =", fixef_include, ", re_include =", re_include, ")")) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 9),
                       axis.title.x = ggplot2::element_text(size = 8),
                       axis.title.y = ggplot2::element_text(size = 8))


      # QQ plot
      qq_vals = qqnorm(res, plot.it = FALSE)
      qq_df = data.frame(Theoretical = qq_vals$x, Sample = qq_vals$y)

      q_sample = quantile(res, probs = c(0.25, 0.75))
      q_theory = qnorm(c(0.25, 0.75))
      qqline_slope = diff(q_sample) / diff(q_theory)
      qqline_intercept = q_sample[1] - qqline_slope * q_theory[1]

      p2 = ggplot2::ggplot(qq_df, ggplot2::aes(x = Theoretical, y = Sample)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(intercept = qqline_intercept, slope = qqline_slope, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste("Normal Q-Q Plot for", plot_title, "Residuals")) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 9),
                       axis.title.x = ggplot2::element_text(size = 8),
                       axis.title.y = ggplot2::element_text(size = 8))

      print(ggpubr::ggarrange(p1, p2, ncol = 2))
    }

    if(model %in% c('downwind_logistic', 'downwind_propensity')){
      df = data.frame(Fitted = fit, Residuals = res)


      #Residuals vs Fitted
      p1 = ggplot2::ggplot(df, ggplot2::aes(x = Fitted, y = Residuals)) +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste("Residuals vs Fitted for", plot_title),
                      y = paste0("Residuals (type=", ifelse(is.null(residual_type), 'deviance', residual_type), ", scaled=", residual_scaled, ")"),
                      x = paste0("Fitted (type=", predict_type, ")") )  +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 9),
                       axis.title.x = ggplot2::element_text(size = 8),
                       axis.title.y = ggplot2::element_text(size = 8))


      # QQ plot
      qq_vals = qqnorm(res, plot.it = FALSE)
      qq_df = data.frame(Theoretical = qq_vals$x, Sample = qq_vals$y)

      q_sample = quantile(res, probs = c(0.25, 0.75))
      q_theory = qnorm(c(0.25, 0.75))
      qqline_slope = diff(q_sample) / diff(q_theory)
      qqline_intercept = q_sample[1] - qqline_slope * q_theory[1]

      p2 = ggplot2::ggplot(qq_df, ggplot2::aes(x = Theoretical, y = Sample)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(intercept = qqline_intercept, slope = qqline_slope, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste("Normal Q-Q Plot for", plot_title, "Residuals")) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 9),
                       axis.title.x = ggplot2::element_text(size = 8),
                       axis.title.y = ggplot2::element_text(size = 8))

      print(ggpubr::ggarrange(p1, p2, ncol = 2))
    }

  }
}


#' @rdname rain_attr-class
#' @export
#'


summary.rain_attr <- function(object, ...) {
  # Extract data name
  data_name <- deparse(object$args$data)

  # Prepare subset expressions with data$var notation
  upwind_subset_expr <- paste(
    deparse(object$args$upwind_subset), " & ", deparse(object$args$positive_subset)
  )
  downwind_subset_expr <- paste(
    deparse(object$args$downwind_subset), " & ", deparse(object$args$positive_subset)
  )

  # Upwind/Downwind LMM
  upwind_fit <- object$all_fitted_models$upwind_lmm_fit
  downwind_fit <- object$all_fitted_models$downwind_lmm_fit

  upwind_lmm_fixef = round(summary(upwind_fit)$coefficients,4)
  upwind_lmm_varcomp = as.data.frame(lme4::VarCorr(upwind_fit))[, c('grp','var1','vcov')]
  colnames(upwind_lmm_varcomp) = c('Groups', 'Name', 'Variance')
  upwind_lmm_varcomp$Name[is.na(upwind_lmm_varcomp$Name)] = ""

  #
  downwind_lmm_fixef = round(summary(downwind_fit)$coefficients,4)
  downwind_lmm_varcomp = as.data.frame(lme4::VarCorr(downwind_fit))[, c('grp','var1','vcov')]
  colnames(downwind_lmm_varcomp) = c('Groups', 'Name', 'Variance')
  downwind_lmm_varcomp$Name[is.na(downwind_lmm_varcomp$Name)] = ""


  # Attribution table
  attr_rows <- c("apo", "apl")
  attr_table <- data.frame(
    Estimate = sapply(attr_rows, function(nm) if(!is.null(object$hatattr[[nm]])) paste0(round(object$hatattr[[nm]]*100,2), "%") else NA),
    row.names = attr_rows
  )

  if(!is.null(object$bootstrap_result)) {
    ci_mat = bootstrap_CI(object$bootstrap_result$hatattr[, attr_rows],level = ifelse(is.null(object$args$bootstrap_option$CI_level), bootstrap_opt()$CI_level, object$args$bootstrap_option$CI_level) )
    attr_table$Bootstrap_CI = apply(ci_mat, 1, function(r) paste0("(", round(r[1] * 100, 4), "%", ", ", round(r[2] * 100, 2), "%", ")"))
    attr_table$Bootstrap_p = round(bootstrap_p_value(object$bootstrap_result$hatattr[, attr_rows]),2)
  }else{
    attr_table$Bootstrap_CI = NA
    attr_table$Bootstrap_p = NA
  }


  if(!is.null(object$permutation_result)) {
    attr_table$Permutation_p = round(permutation_p_value(object$permutation_result$hatattr[, attr_rows],
                                                         ori_est = sapply(attr_rows, FUN = function(nm){
                                                           object$hatattr[[nm]]
                                                         })),2)
  }else{
    attr_table$Permutation_p <- NA
  }

  colnames(attr_table) = c('Estimate',
                           paste0(ifelse(is.null(object$args$bootstrap_option$CI_level), bootstrap_opt()$CI_level, object$args$bootstrap_option$CI_level)*100, "% Bootstrap CI"),
                           'Bootstrap P-Val','Permutation P-Val')


  # SATE table
  sate_rows <- c("sate.mb","sate.ipw","sate.ipw.l","sate.ipw.ma","sate.aipw")
  sate_table <- data.frame(
    Estimate = sapply(sate_rows, function(nm) if(!is.null(object$hatsate[[nm]])) round(object$hatsate[[nm]],4) else NA),
    row.names = sate_rows
  )

  # Bootstrap CI and p-values
  if(!is.null(object$bootstrap_result)) {
    ci_mat_sate <- bootstrap_CI(object$bootstrap_result$hatsate[, sate_rows],
                                level = ifelse(is.null(object$args$bootstrap_option$CI_level),
                                               bootstrap_opt()$CI_level,
                                               object$args$bootstrap_option$CI_level))
    sate_table$Bootstrap_CI <- apply(ci_mat_sate, 1, function(r) paste0("(", round(r[1],4), ", ", round(r[2],4), ")"))
    sate_table$Bootstrap_p <- round(bootstrap_p_value(object$bootstrap_result$hatsate[, sate_rows]),2)
  } else {
    sate_table$Bootstrap_CI <- NA
    sate_table$Bootstrap_p <- NA
  }

  # Permutation p-values
  if(!is.null(object$permutation_result)) {
    sate_table$Permutation_p <- round(permutation_p_value(object$permutation_result$hatsate[, sate_rows],
                                                          ori_est = sapply(sate_rows, function(nm) object$hatsate[[nm]])), 2)
  } else {
    sate_table$Permutation_p <- NA
  }

  colnames(sate_table) = c('Estimate',
                           paste0(ifelse(is.null(object$args$bootstrap_option$CI_level), bootstrap_opt()$CI_level, object$args$bootstrap_option$CI_level)*100, "% Bootstrap CI"),
                           'Bootstrap P-Val','Permutation P-Val')


  summary_list <- list(
    data_name = data_name,

    #upwind
    upwind_subset_expr = upwind_subset_expr,
    upwind_formula = object$args$upwind_lmm_formula,
    upwind_n_obs = nobs(upwind_fit),
    upwind_n_groups = length(unique(lme4::getME(upwind_fit, "flist")[[1]])),
    upwind_summary = summary(upwind_fit),
    upwind_fitted = predict(upwind_fit, re.form = NULL), #including random effects
    upwind_residuals = residuals(upwind_fit, type = 'response', scaled = TRUE),  #including random effects
    upwind_lmm_fixef = upwind_lmm_fixef,
    upwind_lmm_varcomp = upwind_lmm_varcomp,

    #downwind
    downwind_subset_expr = downwind_subset_expr,
    downwind_formula = object$args$downwind_lmm_formula,
    downwind_n_obs= nobs(downwind_fit),
    downwind_n_groups = length(unique(lme4::getME(downwind_fit, "flist")[[1]])),
    downwind_summary = summary(downwind_fit),
    downwind_fitted = predict(downwind_fit, re.form = NULL),  #including random effects
    downwind_residuals = residuals(downwind_fit, type = 'response', scaled = TRUE),  #including random effects
    downwind_lmm_fixef = downwind_lmm_fixef,
    downwind_lmm_varcomp = downwind_lmm_varcomp,

    #
    attr_table = attr_table,
    sate_table = sate_table
  )

  class(summary_list) <- "summary.rain_attr"
  return(summary_list)
}



#' @rdname rain_attr-class
#' @export

print.summary.rain_attr <- function(summary_object, ...) {
  cat("Summary of Two Stage LMM Rainfall Enhancement Analysis\n")
  cat("======================================================================\n\n")

  signif_stars <- function(p) {
    stars <- rep("", length(p))
    stars[!is.na(p) & p < 0.001] <- "***"
    stars[!is.na(p) & p >= 0.001 & p < 0.01] <- "**"
    stars[!is.na(p) & p >= 0.01 & p < 0.05] <- "*"
    stars[!is.na(p) & p >= 0.05 & p < 0.1] <- "."
    return(stars)
  }

  # Attribution Results Table
  cat("Attribution Results (Assuming Log-Rainfall being Modelled):\n")


  attr_tbl = summary_object$attr_table

  # Append stars to all P-Val columns
  pval_cols = grep("P-Val", colnames(attr_tbl))
  for(col in pval_cols){
    stars = signif_stars(as.numeric(attr_tbl[,col]))
    attr_tbl[,col] = paste0(attr_tbl[,col], stars)
  }
  print(attr_tbl)
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")

  # SATE Results Table
  cat("SATE Results:\n")

  sate_tbl <- summary_object$sate_table
  pval_cols <- grep("P-Val", colnames(sate_tbl))
  for(col in pval_cols){
    stars <- signif_stars(as.numeric(sate_tbl[,col]))
    sate_tbl[,col] <- paste0(sate_tbl[,col], stars)
  }
  print(sate_tbl)
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")

  cat("\n======================================================================\n\n")

  # Upwind LMM Results Table
  cat("Upwind (First Stage) LMM:\n")
  cat("Formula:\n")
  print(summary_object$upwind_formula)
  cat('\n')
  cat("Data subset used: ")
  cat(summary_object$data_name, "[", summary_object$upwind_subset_expr, ", ]\n")
  cat(sprintf("Number of observations: %d, Number of groups: %d\n\n",
              summary_object$upwind_n_obs, summary_object$upwind_n_groups))

  # Random effects
  cat("Random effects:\n")
  print(summary_object$upwind_lmm_varcomp, row.names = FALSE)
  cat("\n")

  # Fixed effects
  cat("Fixed effects:\n")
  print(summary_object$upwind_lmm_fixef)
  cat("\n======================================================================\n\n")

  # Downwind LMM Results Table
  cat("Downwind (Second Stage) LMM:\n")
  cat("Formula:\n")
  print(summary_object$downwind_formula)
  cat('\n')
  cat("Data subset used: ")
  cat(summary_object$data_name, "[", summary_object$downwind_subset_expr, ", ]\n")
  cat(sprintf("Number of observations: %d, Number of groups: %d\n\n",
              summary_object$downwind_n_obs, summary_object$downwind_n_groups))

  # Random effects
  cat("Random effects:\n")
  print(summary_object$downwind_lmm_varcomp, row.names = FALSE)
  cat("\n")

  # Fixed effects
  cat("Fixed effects:\n")
  print(summary_object$downwind_lmm_fixef)

  invisible(summary_object)

}
