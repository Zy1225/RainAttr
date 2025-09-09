
#TODO: Add the following types of EDA

# Table summarizing number of observations/unique groups of Upwind/Downwind/Target/Control vs. Positive/Zero (DONE)

# Histogram of cluster sizes for Upwind+Positive / Downwind+Positive

# QQplots of raw/log rainfall for Upwind+Positive / Downwind+Positive

# 2-dimensional heatmap plot for raw/lograinfall (including those with zeros), with row = gauge, column = day

# Time series plot (separate for each year) with each line representing mean (with and without zeros) raw/log rainfall of Upwind/Target/Control

# Time series plot (separate for each year) with each line representing the OBSERVED raw/log rainfall of each gauge, with optional filtering of gauges

# A spatial plot of map (separate for each year) with dots colored based on the mean raw/log rainfall (averaged within each year, with and without zeros) of each gauge


eda = function(eda_type,
               data, rain_col_name, day_column_name, year_column_name, use_raw, filter_gauge = NULL,
               upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset, positive_subset){

  original_args = as.list(match.call())[-1]

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
  downwind_target_expr = rlang::enquo(downwind_target_subset)
  target = rlang::eval_tidy(downwind_target_expr, data = data)

  #Binary indicator of length N_downwind_positive, indicating whether or not each downwind positive observation is control observation
  downwind_control_expr = rlang::enquo(downwind_control_subset)
  control = rlang::eval_tidy(downwind_control_expr, data = data)

  if(eda_type == 'num_obs_days'){
    output_obs = rbind(
      c(sum(upwind & positive), sum(upwind & !positive)),
      c(sum(downwind & positive), sum(downwind & !positive)),
      c(sum(target & positive), sum(target & !positive)),
      c(sum(control & positive), sum(control & !positive))
    )

    output_days = rbind(
      c(length(unique(data[upwind & positive, day_column_name])) , length(unique(data[upwind & !positive, day_column_name]))),
      c(length(unique(data[downwind & positive, day_column_name])) , length(unique(data[downwind & !positive, day_column_name]))),
      c(length(unique(data[target & positive, day_column_name])) , length(unique(data[target & !positive, day_column_name]))),
      c(length(unique(data[control & positive, day_column_name])) , length(unique(data[control & !positive, day_column_name])))
    )

    rownames(output_obs) = rownames(output_days) = c(original_args$upwind_subset, original_args$downwind_subset, original_args$downwind_target_subset, original_args$downwind_control_subset)
    colnames(output_obs) = colnames(output_days) = c(original_args$positive_subset, paste0("!(", deparse(original_args$positive_subset), ")"))

    return(list(
      num_obs = output_obs,
      num_unique_days = output_days
    ))
  }



  if(eda_type == 'hist_day_group_sizes'){
    day_values = data[[day_column_name]]

    # Upwind & Positive group sizes
    upwind_days = day_values[upwind & positive]
    upwind_sizes = as.numeric(table(upwind_days))

    # Downwind (Target + Control) & Positive group sizes
    downwind_days = day_values[downwind & positive]
    downwind_sizes = as.numeric(table(downwind_days))

    upwind_title <- paste0("Histogram of Day Group Sizes \n(",
                           deparse(original_args$upwind_subset),
                           " & ",
                           deparse(original_args$positive_subset), ")")

    downwind_title <- paste0("Histogram of Day Group Sizes \n(",
                             deparse(original_args$downwind_subset),
                             " & ",
                             deparse(original_args$positive_subset), ")")

    # Create histograms
    p1 <- ggplot2::ggplot(data.frame(size = upwind_sizes), ggplot2::aes(x = size)) +
      ggplot2::geom_histogram(binwidth = 1, fill = "#1f77b4", color = "black") +
      ggplot2::labs(title = upwind_title, x = "Group Size", y = "Frequency") +
      ggplot2::theme_bw()

    p2 <- ggplot2::ggplot(data.frame(size = downwind_sizes), ggplot2::aes(x = size)) +
      ggplot2::geom_histogram(binwidth = 1, fill = "#ff7f0e", color = "black") +
      ggplot2::labs(title = downwind_title, x = "Group Size", y = "Frequency") +
      ggplot2::theme_bw()

    print(ggpubr::ggarrange(plotlist = list(p1,p2)))

    return(list(
      upwind_positive_hist = p1,
      downwind_positive_hist = p2
    ))
  }
}

qwe =eda(eda_type = 'hist_day_group_sizes',
    data = oman,


    #Specify the column in input data that contains the raw rainfall
    rain_col_name = 'Rain.Gauge.Measurement',
    day_column_name = 'TrialDay',

    #Logical expression identifying subset of observations to which upwind_lmm_formula is fitted
    upwind_subset = Gauge.Day.Type == 'Upwind',

    #Logical expression identifying subset of observations to which downwind_lmm_formula is fitted
    downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),

    #Logical expression identifying subset of observations of those that are 'Target'
    downwind_target_subset = Gauge.Day.Type == 'Target',

    #Logical expression identifying subset of observations of those that are 'Control'
    downwind_control_subset = Gauge.Day.Type == 'Control',

    #Logical expression identifying subset of observations with rainfall event
    positive_subset = Rain.Gauge.Measurement > 0,
)
