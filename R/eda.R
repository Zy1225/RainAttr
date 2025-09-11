
#TODO: Add the following types of EDA

# Table summarizing number of observations/unique groups of Upwind/Downwind/Target/Control vs. Positive/Zero (DONE)
# Maybe add separate tables for each year? (DONE)

# Histogram of cluster sizes for Upwind+Positive / Downwind+Positive (DONE)
# Maybe add separate plots for each year?
# - maybe not since we are not really interested in each year's cluster size distribution
# - Also, we only really want to look at the OVERALL cluster size distribution, to decide what bootstrap to be used, since we are fitting model to all years but not separately to each year
# - Also, we could generate separate plots for each year by manually hacking eda(), by setting positive_subset = Rain.Gauge.Measurement > 0 & Year == XXX

# QQplots of raw/log rainfall for Upwind+Positive / Downwind+Positive (DONE)
# Maybe add separate plots for each year?
# - Maybe not since we are fitting models to ALL years but not separately to each year
# - Also, we could generate separate plots for each year by manually hacking eda(), by setting positive_subset = Rain.Gauge.Measurement > 0 & Year == XXX

# 2-dimensional heatmap plot for raw/lograinfall (including those with zeros), with row = gauge, column = day
# Maybe not since this is not really useful due to not having spatial information

# Time series plot (separate for each year) with each line representing mean (with and without zeros) raw/log rainfall of Upwind/Target/Control (DONE)
#Maybe just consider mean without zeros, since we are modelling only non-zeros as well as there are TOO MANY ZEROS

# Time series plot (separate for each year) with each line representing the OBSERVED raw/log rainfall of each gauge, with optional filtering of gauges
#Maybe need to add a restriction to only plot a small number of gauges

# A static spatial plot of map (separate for each year) with dots colored based on the mean raw/log rainfall (averaged within each year, with and without zeros) of each gauge

#A dynamic spatial plot of map (focused on each year) with dots colored based on the daily raw/log rainfall of each gauge
#using gganimate and transition_time(day)

#Think if we want to add smoothed (static/dynamic) spatial map

#TODO: For documentation, remember to add a note that log(0) are plotted at the bottommost of the plot
#TODO: Add another version of ts_by_gauge_interactive to enable ts_focus_gauge
#TODO: Check if ts_by_gauge_interactive's plotly viewer is working correctly when we have A LARGE number of gauges
#TODO: Add another version of map_dynamic when input_sf is NULL, and also another version when chosen_year=NULL
eda = function(eda_type,
               data, rain_col_name, day_column_name, year_column_name, use_raw,
               gauge_id_column_name, ts_focus_gauge = NULL,
               longlat_column_names, long_lim, lat_lim, input_sf = NULL,
               chosen_year,
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

  if(eda_type == 'num_obs_days_by_year'){
    final_output =
      lapply(
        unique(data[,year_column_name]), FUN = function(year_val){
          year_ind = (data[,year_column_name] == year_val)

          output_obs = rbind(
            c(sum(year_ind & upwind & positive), sum(year_ind & upwind & !positive)),
            c(sum(year_ind & downwind & positive), sum(year_ind & downwind & !positive)),
            c(sum(year_ind & target & positive), sum(year_ind & target & !positive)),
            c(sum(year_ind & control & positive), sum(year_ind & control & !positive))
          )

          output_days = rbind(
            c(length(unique(data[year_ind & upwind & positive, day_column_name])) , length(unique(data[year_ind & upwind & !positive, day_column_name]))),
            c(length(unique(data[year_ind & downwind & positive, day_column_name])) , length(unique(data[year_ind & downwind & !positive, day_column_name]))),
            c(length(unique(data[year_ind & target & positive, day_column_name])) , length(unique(data[year_ind & target & !positive, day_column_name]))),
            c(length(unique(data[year_ind & control & positive, day_column_name])) , length(unique(data[year_ind & control & !positive, day_column_name])))
          )

          rownames(output_obs) = rownames(output_days) = c(original_args$upwind_subset, original_args$downwind_subset, original_args$downwind_target_subset, original_args$downwind_control_subset)
          colnames(output_obs) = colnames(output_days) = c(original_args$positive_subset, paste0("!(", deparse(original_args$positive_subset), ")"))

          return(list(
            num_obs = output_obs,
            num_unique_days = output_days
          ))
        }
      )

    names(final_output) = unique(data[,year_column_name])

    return(final_output)
  }



  if(eda_type == 'hist_day_group_sizes'){
    day_values = data[,day_column_name]

    # Upwind & Positive group sizes
    upwind_days = day_values[upwind & positive]
    upwind_sizes = as.numeric(table(upwind_days))

    # Downwind (Target + Control) & Positive group sizes
    downwind_days = day_values[downwind & positive]
    downwind_sizes = as.numeric(table(downwind_days))

    upwind_title <- paste0("Histogram for Day Group Sizes of\n(",
                           deparse(original_args$upwind_subset),
                           ") & (",
                           deparse(original_args$positive_subset), ")")

    downwind_title <- paste0("Histogram for Day Group Sizes of\n(",
                             deparse(original_args$downwind_subset),
                             ") & (",
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

  if(eda_type == 'qq_rain'){

    if(use_raw){
      rain_values = data[, rain_col_name]
      rain_label = rain_col_name
    }else{
      rain_values = log(data[, rain_col_name])
      rain_label = paste0("log(", rain_col_name, ")")
    }

    # QQplots of raw/log rainfall for Upwind+Positive / Downwind+Positive

    qq_vals_upwind = qqnorm(rain_values[upwind & positive], plot.it = FALSE)
    qq_df_upwind = data.frame(Theoretical = qq_vals_upwind$x, Sample = qq_vals_upwind$y)

    q_sample_upwind = quantile(rain_values[upwind & positive], probs = c(0.25, 0.75))
    q_theory_upwind = qnorm(c(0.25, 0.75))
    qqline_slope_upwind = diff(q_sample_upwind) / diff(q_theory_upwind)
    qqline_intercept_upwind = q_sample_upwind[1] - qqline_slope_upwind * q_theory_upwind[1]


    p1 = ggplot2::ggplot(qq_df_upwind, ggplot2::aes(x = Theoretical, y = Sample)) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(intercept = qqline_intercept_upwind, slope = qqline_slope_upwind, linetype = "dashed", color = "red") +
      ggplot2::labs(title = paste0("Normal Q-Q Plot for ", rain_label, " of\n (",
                                  deparse(original_args$upwind_subset),
                                  ") & (",
                                  deparse(original_args$positive_subset), ")")) +
      ggplot2::theme_bw()


    qq_vals_downwind = qqnorm(rain_values[downwind & positive], plot.it = FALSE)
    qq_df_downwind = data.frame(Theoretical = qq_vals_downwind$x, Sample = qq_vals_downwind$y,
                                Group = factor(ifelse(target[downwind & positive], "Target", "Control"))
    )

    q_sample_downwind = quantile(rain_values[downwind & positive], probs = c(0.25, 0.75))
    q_theory_downwind = qnorm(c(0.25, 0.75))
    qqline_slope_downwind = diff(q_sample_downwind) / diff(q_theory_downwind)
    qqline_intercept_downwind = q_sample_downwind[1] - qqline_slope_downwind * q_theory_downwind[1]

    legend_title = as.character(original_args$downwind_target_subset)[2]

    p2 = ggplot2::ggplot(qq_df_downwind, ggplot2::aes(x = Theoretical, y = Sample, color = Group, shape = Group)) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(intercept = qqline_intercept_downwind, slope = qqline_slope_downwind, linetype = "dashed", color = "red") +
      ggplot2::labs(title = paste0("Normal Q-Q Plot for ", rain_label, "of \n(",
                                  deparse(original_args$downwind_subset),
                                  ") & (",
                                  deparse(original_args$positive_subset), ")"),
                    color = legend_title, shape = legend_title) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") +
      ggplot2::scale_color_manual(values = c("Target" = "#1f77b4", "Control" = "#ff7f0e")) +
      ggplot2::scale_shape_manual(values = c("Target" = 16, "Control" = 17))

    print(ggpubr::ggarrange(plotlist = list(p1,p2)))

    return(list(
      upwind_positive_qq = p1,
      downwind_positive_qq = p2
    ))
  }

  if(eda_type == 'ts_by_type'){

    if(!use_raw){
      data[,rain_col_name] = log(data[,rain_col_name])
      rain_label = paste0("log(", rain_col_name, ")")
    }else{
      rain_label = rain_col_name
    }



    upwind_df = data[upwind & positive, c(day_column_name, year_column_name, rain_col_name)]
    upwind_df_avg = aggregate(as.formula(
      paste(rain_col_name, "~", paste(c(day_column_name, year_column_name), collapse = " + "))
    ), data = upwind_df, FUN = mean, na.rm = TRUE)
    upwind_df_avg = cbind(upwind_df_avg, type = 'upwind')

    target_df = data[target & positive, c(day_column_name, year_column_name, rain_col_name)]
    target_df_avg = aggregate(as.formula(
      paste(rain_col_name, "~", paste(c(day_column_name, year_column_name), collapse = " + "))
    ), data = target_df, FUN = mean, na.rm = TRUE)
    target_df_avg = cbind(target_df_avg, type = 'target')

    control_df = data[control & positive, c(day_column_name, year_column_name, rain_col_name)]
    control_df_avg = aggregate(as.formula(
      paste(rain_col_name, "~", paste(c(day_column_name, year_column_name), collapse = " + "))
    ), data = control_df, FUN = mean, na.rm = TRUE)
    control_df_avg = cbind(control_df_avg, type = 'control')


    final_df_avg = rbind(upwind_df_avg, target_df_avg, control_df_avg)
    final_df_avg$type = factor(final_df_avg$type, levels = c('target','control','upwind'))
    #browser()

    output_plot = ggplot2::ggplot(final_df_avg, ggplot2::aes(x = .data[[day_column_name]], y = .data[[rain_col_name]], group = type,
                                                             color = type)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() + ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]]), scales = 'free') +
      ggplot2::labs(title = paste0('Time Series Plots of ', rain_label, ' Averaged Across\nDifferent Types of Observations Satisfying ',
                                   deparse(original_args$positive_subset)),
                    y = rain_label,
                    color = NULL) +
      ggplot2::scale_color_manual(values = c("target" = "#1f77b4", "control" = "#ff7f0e", "upwind" = "#9467bd"),
                                  labels = c(deparse(original_args$downwind_target_subset),
                                             deparse(original_args$downwind_control_subset),
                                             deparse(original_args$upwind_subset))) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")


    print(output_plot)

    return(output_plot)

  }

  if(eda_type == 'ts_by_gauge'){
    if(!use_raw){
      data[,rain_col_name] = log(data[,rain_col_name])
      rain_label = paste0("log(", rain_col_name, ")")
    }else{
      rain_label = rain_col_name
    }

    if(!is.null(ts_focus_gauge)){
      data[,gauge_id_column_name] = factor(data[,gauge_id_column_name])
      output_plot = ggplot2::ggplot(data, ggplot2::aes(x = .data[[day_column_name]], y = .data[[rain_col_name]], group = .data[[gauge_id_column_name]])) +
        ggplot2::geom_line(data = data[!data[,gauge_id_column_name] %in% ts_focus_gauge,] ,
                           ggplot2::aes(color = "Other"), alpha = 0.8) +
        ggplot2::geom_point(data = data[!data[,gauge_id_column_name] %in% ts_focus_gauge,] ,
                            ggplot2::aes(color = "Other"), alpha = 0.8) +
        ggplot2::geom_line(data = data[data[,gauge_id_column_name] %in% ts_focus_gauge,] ,
                           ggplot2::aes(color = .data[[gauge_id_column_name]]), alpha = 0.8, linewidth = 1) +
        ggplot2::geom_point(data = data[data[,gauge_id_column_name] %in% ts_focus_gauge,] ,
                            ggplot2::aes(color = .data[[gauge_id_column_name]]), alpha = 0.8) +
        ggplot2::scale_color_manual(
          values = c("Other" = "grey80",
                     setNames(viridis::viridis(length(ts_focus_gauge), option = "D"),
                              ts_focus_gauge))
        ) +
        ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]]), scales = 'free') +
        ggplot2::labs(title = paste0('Time Series Plots of ', rain_label, ' of All Gauges'),
                      y = rain_label,
                      color = gauge_id_column_name) +
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = "bottom")
    }else{
      output_plot = ggplot2::ggplot(data, ggplot2::aes(x = .data[[day_column_name]], y = .data[[rain_col_name]], group = .data[[gauge_id_column_name]])) +
        ggplot2::geom_line(alpha = 0.8, color = "grey80") +
        ggplot2::geom_point(color = "grey80") +
        ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]]), scales = 'free') +
        ggplot2::labs(title = paste0('Time Series Plots of ', rain_label, ' of All Gauges'),
                      y = rain_label,
                      color = gauge_id_column_name) +
        ggplot2::theme_bw()
    }


    print(output_plot)

    return(output_plot)


  }

  if(eda_type == 'ts_by_gauge_interactive'){
    if(!use_raw){
      data[,rain_col_name] = log(data[,rain_col_name])
      rain_label = paste0("log(", rain_col_name, ")")
    }else{
      rain_label = rain_col_name
    }


    data[,gauge_id_column_name] <- factor(data[,gauge_id_column_name])

    output_plot = ggplot2::ggplot(data, ggplot2::aes(x = .data[[day_column_name]],
                                                     y = .data[[rain_col_name]],
                                                     group = .data[[gauge_id_column_name]],
                                                     # color = .data[[gauge_id_column_name]],
                                                     text = paste("Gauge:", .data[[gauge_id_column_name]]))) +
      ggplot2::geom_line(alpha = 0.8, color = "grey80") +
      ggplot2::geom_point(color = "grey80") +
      ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]]), scales = "free") +
      ggplot2::labs(title = paste0('Time Series Plots of ', rain_label, ' of All Gauges (Hover over points to see Gauge ID)'),
                    y = rain_label) +
      ggplot2::theme_bw()

    plotly_output = plotly::ggplotly(output_plot, tooltip = "text")

    print(plotly_output)
    print('hi')
    return(plotly_output)
  }

  if(eda_type == "map_static"){
    if(!use_raw){
      data[,rain_col_name] = log(data[,rain_col_name])
      rain_label = paste0("log(", rain_col_name, ")")
    }else{
      rain_label = rain_col_name
    }

    positive_df = data[positive, c(longlat_column_names, year_column_name, rain_col_name)]
    positive_df_avg = aggregate(as.formula(
      paste(rain_col_name, "~", paste(c(longlat_column_names, year_column_name), collapse = " + "))
    ), data = positive_df, FUN = mean, na.rm = TRUE)


    #the following need either 'maps' package or 'rnaturalearthhires' for creating input_sf
    if(!is.null(input_sf)){
      points_sf <- sf::st_as_sf(
        positive_df_avg,
        coords = c(longlat_column_names[1], longlat_column_names[2]),
        crs = 4326
      )

      output_plot = ggplot2::ggplot() +
        ggplot2::geom_sf(data = input_sf, fill = "white", color = "black") +
        ggplot2::geom_sf(data = points_sf, ggplot2::aes(color = .data[[rain_col_name]])) +
        ggplot2::coord_sf(xlim = long_lim, ylim = lat_lim, expand = FALSE) +
        ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]])) +
        ggplot2::labs(title = paste0('Spatial Plots of ', rain_label, ' Averaged Across\n All Days Satisfying ',
                                     deparse(original_args$positive_subset)),
                      x = 'Longitude', y = 'Latitude', color = rain_label) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = 'bottom')

    }else{
      output_plot = ggplot2::ggplot(positive_df_avg, ggplot2::aes(x = .data[[longlat_column_names[1]]],
                                                                  y = .data[[longlat_column_names[2]]],
                                                                  color = .data[[rain_col_name]])) +
        ggplot2::borders()+
        ggplot2::xlim(long_lim) + ggplot2::ylim(lat_lim) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]]), scales = "free") +
        ggplot2::labs(title = paste0('Spatial Plots of ', rain_label, ' Averaged Across\n All Days Satisfying ',
                                     deparse(original_args$positive_subset)),
                      x = 'Longitude', y = 'Latitude', color = rain_label) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = 'bottom')
    }

    print(output_plot)
    return(output_plot)


  }

  if(eda_type == 'map_dynamic'){

    if(!use_raw){
      data[,rain_col_name] = log(data[,rain_col_name])
      rain_label = paste0("log(", rain_col_name, ")")
    }else{
      rain_label = rain_col_name
    }

    positive_df = data[positive, c(longlat_column_names, day_column_name, year_column_name, rain_col_name)]

    points_sf <- sf::st_as_sf(
      positive_df,
      coords = c(longlat_column_names[1], longlat_column_names[2]),
      crs = 4326
    )


    points_sf_year <- points_sf[points_sf[[year_column_name]] == chosen_year, ]

    output_plot = ggplot2::ggplot() +
      ggplot2::geom_sf(data = input_sf, fill = "white", color = "black", linewidth = 0.3) +
      ggplot2::geom_sf(
        data = points_sf_year,
        ggplot2::aes(color = .data[[rain_col_name]]),
        size = 2
      ) +
      ggplot2::coord_sf(xlim = long_lim, ylim = lat_lim, expand = FALSE) +
      ggplot2::labs(
        title = paste0(
          "Spatial Plots of ", rain_label, " for Year ", chosen_year,
          "\nDay: {closest_state}"
        ),
        x = "Longitude",
        y = "Latitude",
        color = rain_label
      ) +
      ggplot2::scale_color_viridis_c(name = rain_label) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom") +
      gganimate::transition_states(
        states = .data[[day_column_name]],
        state_length = 1,
        transition_length = 1
      ) +
      gganimate::ease_aes("linear")

    gganimate::animate(output_plot)
    return(gganimate::animate(output_plot))
  }

}

qwe =eda(eda_type = 'map_dynamic',
         # 'num_obs_days'
         # 'num_obs_days_by_year'
         # 'hist_day_group_sizes'
         # 'qq_rain'
    # data = oman,
    data = oman,


    #Specify the column in input data that contains the raw rainfall
    rain_col_name = 'Rain.Gauge.Measurement',
    day_column_name = 'TrialDay',
    year_column_name = 'Year',

    use_raw = F,

    gauge_id_column_name = 'Gauge.ID', ts_focus_gauge = c(1,3,212,213),
    longlat_column_names = c("Gauge.Longitude", "Gauge.Latitude"),
    long_lim = c(55,60),
    lat_lim = c(20,25),
    input_sf = rnaturalearth::ne_countries(scale = "large", country = "Oman", returnclass = "sf"),

    chosen_year = 2015,

    #Logical expression identifying subset of observations to which upwind_lmm_formula is fitted
    upwind_subset = Gauge.Day.Type == 'Upwind',

    #Logical expression identifying subset of observations to which downwind_lmm_formula is fitted
    downwind_subset = Gauge.Day.Type  %in% c('Target','Control'),

    #Logical expression identifying subset of observations of those that are 'Target'
    downwind_target_subset = Gauge.Day.Type == 'Target',

    #Logical expression identifying subset of observations of those that are 'Control'
    downwind_control_subset = Gauge.Day.Type == 'Control',

    #Logical expression identifying subset of observations with rainfall event
    positive_subset = Rain.Gauge.Measurement > 0
)
# qwe
#do.call(sum,lapply(qwe, function(x){x$num_obs}))
