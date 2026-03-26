#' Exploratory Data Analysis for Rainfall Enhancement Trial Data
#'
#' Performs various exploratory data analyses (EDA) on rainfall enhancement trial data,
#' including summaries of observation types and positive rainfall events, time series plots, Q-Q plots, and static or animated spatial maps.
#'
#' @param eda_type A character string specifying the type of EDA to perform. Must be one of
#'   \code{"num_obs_days"}, \code{"num_obs_days_by_year"}, \code{"hist_day_group_sizes"},
#'   \code{"qq_rain"}, \code{"ts_by_type"}, \code{"ts_by_gauge"},
#'   \code{"ts_by_gauge_interactive"}, \code{"map_static"}, \code{"map_dynamic"}.
#' @param data A data frame containing the rainfall enhancement trial data.
#' @param rain_col_name A character string that refers to the column name of the raw scale rainfall in \code{data}.
#' @param day_column_name A character string that refers to the column name of day identifiers in \code{data}.
#' @param year_column_name A character string that refers to the column name of year identifiers in \code{data}.
#' @param use_raw Logical. If TRUE, EDA is performed on raw scale rainfall. if FALSE, EDA is performed on log-transformed rainfall.
#' @param gauge_id_column_name A character string that refers to the column name of gauge identifiers in \code{data}.
#' @param ts_focus_gauge An optional vector of gauge identifiers to highlight in time series plots. If \code{ts_focus_gauge} is not supplied, no highlighting is done and all gauges are plotted with the same color.
#' @param longlat_column_names A character vector of length 2 specifying the column names of longitude and latitude in \code{data}, for plotting spatial plots.
#' @param long_lim A numeric vector of length 2 specifying longitude limits for spatial plots.
#' @param lat_lim A numeric vector of length 2 specifying latitude limits for spatial plots.
#' @param input_sf An optional \code{\link[sf:st_as_sf]{sf}} object for polygon layers (e.g., country borders) in spatial plots. For example, \code{rnaturalearth::ne_countries(scale = "large", country = "Oman", returnclass = "sf")}.
#' If \code{input_sf} is not supplied, then a default map is drawn using \code{ggplot2::borders()}.
#' @param ionizer_location_df An optional data frame containing the locations of ionizers, to be used in the plotting of static and animated maps.
#' If \code{ionizer_location_df} is not supplied, no ionizers are plotted.
#' The data frame must contain columns specified by \code{ionizer_id_column_name} and \code{ionizer_longlat_column_names}. It should also include additional columns whose names correspond to the years in \code{unique(data[[year_column_name]])}, with binary values indicating whether each ionizer has been deployed in each year.
#' @param ionizer_id_column_name An optional character string that refers to the column name of ionizer identifiers in \code{ionizer_location_df}. This must be supplied if \code{ionizer_location_df} is supplied.
#' @param ionizer_longlat_column_names An optional character vector of length 2 specifying the column names of longitude and latitude in \code{ionizer_location_df}, for plotting ionizers. This must be supplied if \code{ionizer_location_df} is supplied.
#' @param elev_contour An optional logical. If \code{TRUE}, elevation contour lines are added to the spatial plots. If \code{FALSE} (default), elevation contour lines are not included.
#' @param elev_resolution An optional integer between 1 and 14, specifying the resolution of the elevation data obtained from Amazon Web Services Terrain Tiles via \code{\link[elevatr]{get_elev_raster}}. Defaults to 2. Higher values indicate higher resolution; see the \code{z} argument in \code{\link[elevatr]{get_elev_raster}} for more details.
#' @param focus_year An optional vector specifying years to filter for animated maps. If \code{focus_year} is not supplied, all years are included.
#' @param fps An optional numeric specifying frames per second for animated maps. Default is \code{10}.
#' @param animate_filename An optional character string specifying the file name used to save the animated map as a GIF. The file name should end with \code{".gif"}. If \code{animate_filename} is not supplied, no GIF is saved on disk.
#' @param upwind_subset A logical expression used to extract the relevant subset of observations from \code{data} to be used in the upwind (first stage) LMM fitting. For example, \code{Gauge.Day.Type == "Upwind"}.
#' @param downwind_subset A logical expression used to extract the relevant subset of observations from \code{data} to be used in the downwind (second stage) LMM fitting. For example, \code{Gauge.Day.Type \%in\% c("Target","Control")}.
#' @param downwind_target_subset A logical expression used to extract the relevant subset of downwind (second stage) observations from \code{data} that were exposed to treatment (operating ionizers). For example, \code{Gauge.Day.Type == "Target"}.
#' @param downwind_control_subset A logical expression used to extract the relevant subset of downwind (second stage) observations from \code{data} that were not exposed to treatment (operating ionizers). For example, \code{Gauge.Day.Type == "Control"}.
#' @param positive_subset A logical expression used to extract the relevant subset of observations from \code{data} with positive rainfall. For example, \code{Rain.Gauge.Measurement > 0}.
#'
#' @details
#' Each \code{eda_type} has its own behavior and relevant arguments provided in parentheses:
#'
#' \describe{
#'   \item{`num_obs_days` (\code{data}, \code{day_column_name}, \code{upwind_subset}, \code{downwind_subset}, \code{downwind_target_subset}, \code{downwind_control_subset}, \code{positive_subset})}{
#'     Computes contingency tables of the number of observations and number of unique days for each subset
#'     (\code{upwind_subset}, \code{downwind_subset}, \code{downwind_target_subset}, \code{downwind_control_subset}), with columns corresponding to positive vs zero rainfall events.
#'   }
#'   \item{`num_obs_days_by_year` (\code{data}, \code{day_column_name}, \code{year_column_name}, \code{upwind_subset}, \code{downwind_subset}, \code{downwind_target_subset}, \code{downwind_control_subset}, \code{positive_subset})}{
#'     Same as \code{num_obs_days}, but computed separately for each year.
#'   }
#'   \item{`hist_day_group_sizes` (\code{data}, \code{day_column_name}, \code{upwind_subset}, \code{downwind_subset}, \code{positive_subset})}{
#'     Plots histograms of group sizes for days in \code{upwind_subset} and \code{downwind_subset} with positive rainfall.
#'     The group size for a given day and type (\code{upwind_subset}, \code{downwind_subset}) is defined as the number of gauges satisfying that type and having positive rainfall on that day.
#'   }
#'   \item{`qq_rain` (\code{data}, \code{rain_col_name}, \code{use_raw}, \code{upwind_subset}, \code{downwind_subset}, \code{downwind_target_subset}, \code{downwind_control_subset}, \code{positive_subset})}{
#'     Produces Normal Q-Q plots for rainfall values (raw or log-transformed) for all observations satisfying \code{upwind_subset & positive_subset}, as well as those observations satisfying \code{downwind_subset & positive_subset}.
#'   }
#'   \item{`ts_by_type` (\code{data}, \code{rain_col_name}, \code{day_column_name}, \code{year_column_name}, \code{use_raw}, \code{upwind_subset}, \code{downwind_target_subset}, \code{downwind_control_subset}, \code{positive_subset})}{
#'     Plots daily average rainfall (raw or log-transformed) by subset (\code{upwind_subset}, \code{downwind_target_subset}, \code{downwind_control_subset}), facetted by year.
#'     Averaging is performed only over observations with positive rainfall within each day.
#'   }
#'   \item{`ts_by_gauge` (\code{data}, \code{rain_col_name}, \code{day_column_name}, \code{year_column_name}, \code{use_raw}, \code{gauge_id_column_name}, \code{ts_focus_gauge})}{
#'     Plots daily rainfall (raw or log-transformed) time series for each gauge, optionally highlighting a subset of gauges, with faceting by year.
#'     When plotting log-transformed rainfall using \code{\link[ggplot2:ggplot]{ggplot2}}, observations with zero rainfall are represented as a point at the bottommost of the plot.
#'   }
#'   \item{`ts_by_gauge_interactive` (\code{data}, \code{rain_col_name}, \code{day_column_name}, \code{year_column_name}, \code{use_raw}, \code{gauge_id_column_name}, \code{ts_focus_gauge})}{
#'     Interactive version of \code{ts_by_gauge} using \code{\link[plotly:plot_ly]{plotly}}. Points show gauge identifiers on hover.
#'     When plotting log-transformed rainfall (\code{use_raw = FALSE}), zero rainfall values result in \code{-Inf}.
#'     Unlike \code{\link[ggplot2:ggplot]{ggplot2}}, \code{\link[plotly:plot_ly]{plotly}} will omit these points, breaking the lines. Therefore, users should be aware that lines may appear broken for days with zero rainfall when plotting log-transformed rainfall.
#'   }
#'   \item{`map_static` (\code{data}, \code{rain_col_name}, \code{day_column_name}, \code{year_column_name}, \code{use_raw},  \code{longlat_column_names}, \code{long_lim}, \code{lat_lim}, \code{input_sf}, \code{ionizer_location_df}, \code{ionizer_id_column_name}, \code{ionizer_longlat_column_names}, \code{elev_contour}, \code{elev_resolution}, \code{positive_subset})}{
#'     Produces a static spatial map of annual average rainfall (raw or log-transformed), optionally overlaying an \code{\link[sf:st_as_sf]{sf}} polygon layer supplied via \code{input_sf} and adding elevation contour lines if \code{elev_contour = TRUE} as well as plotting ionizers supplied via \code{ionizer_location_df}, with faceting by year.
#'     Averaging is performed only over days with positive rainfall for each gauge. Requires the \pkg{maps} package to draw map borders when \code{input_sf} is not supplied. Also requires the \pkg{elevatr} and \pkg{raster} package to plot elevation contour lines when \code{elev_contour = TRUE}.
#'   }
#'   \item{`map_dynamic` (\code{data}, \code{rain_col_name}, \code{day_column_name}, \code{year_column_name},  \code{use_raw}, \code{longlat_column_names}, \code{long_lim}, \code{lat_lim}, \code{input_sf}, \code{ionizer_location_df}, \code{ionizer_id_column_name}, \code{ionizer_longlat_column_names}, \code{elev_contour}, \code{elev_resolution}, \code{focus_year}, \code{fps}, \code{animate_filename})}{
#'     Produces an animated map showing rainfall (raw or log-transformed) for each day, optionally filtered by year using the argument \code{focus_year} and overlaying an \code{\link[sf:st_as_sf]{sf}} polygon layer supplied via \code{input_sf} as well as adding elevation contour lines if \code{elev_contour = TRUE}, and plotting ionizers supplied via \code{ionizer_location_df}.
#'     The animation frames are displayed in the order of \code{data[,day_column_name]}. Users should ensure that \code{data[,day_column_name]} contains values that can be meaningfully ordered (e.g., numeric or Date), rather than nominal/factor values, so the animation reflects the correct temporal progression.
#'     Requires the \pkg{maps} package to draw map borders when \code{input_sf} is not supplied and the \pkg{gifski} package for rendering animated map. Also requires the \pkg{elevatr} and \pkg{raster} package to plot elevation contour lines when \code{elev_contour = TRUE}. When \code{animate_filename} is supplied, the resulting gif_image of the animated map will be saved to the location specified by \code{animate_filename}.
#'   }
#' }
#'
#' @return Depending on \code{eda_type}:
#'   - Data summaries as lists (\code{"num_obs_days"}, \code{"num_obs_days_by_year"})
#'   - \code{\link[ggplot2]{ggplot}} objects (\code{"hist_day_group_sizes"}, \code{"qq_rain"}, \code{"ts_by_type"},
#'     \code{"ts_by_gauge"}, \code{"map_static"})
#'   - Interactive \code{\link[plotly:plot_ly]{plotly}} object (\code{"ts_by_gauge_interactive"})
#'   - A list consisting of an animated \code{\link[gganimate:gganimate-package]{gganim}} object and its associated gif_image (\code{"map_dynamic"})
#'


#TODO: Think if we want to add smoothed (static/dynamic) spatial map

eda = function(eda_type,
               data, rain_col_name, day_column_name, year_column_name, use_raw,
               gauge_id_column_name, ts_focus_gauge = NULL,
               longlat_column_names, long_lim, lat_lim, input_sf = NULL,
               ionizer_location_df = NULL, ionizer_id_column_name = NULL, ionizer_longlat_column_names = NULL,
               elev_contour = TRUE, elev_resolution = 2,
               focus_year = NULL, fps = 10, animate_filename = NULL,
               upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset, positive_subset){

  original_args = as.list(match.call())[-1]


  allowed_eda_type = c("num_obs_days", "num_obs_days_by_year", "hist_day_group_sizes", "qq_rain", "ts_by_type", "ts_by_gauge", "ts_by_gauge_interactive", "map_static", "map_dynamic")
  if(!eda_type %in% allowed_eda_type){
    stop(paste("Invalid eda_type. Choose one of:", paste(allowed_eda_type, collapse = ", ")))
  }

  if (eda_type == "num_obs_days") {
    if (missing(data) || missing(day_column_name) ||
        missing(upwind_subset) || missing(downwind_subset) ||
        missing(downwind_target_subset) || missing(downwind_control_subset) ||
        missing(positive_subset)) {
      stop("When eda_type = 'num_obs_days', you must supply: data, day_column_name, upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset, positive_subset.")
    }
  }

  if (eda_type == "num_obs_days_by_year") {
    if (missing(data) || missing(day_column_name) || missing(year_column_name) ||
        missing(upwind_subset) || missing(downwind_subset) ||
        missing(downwind_target_subset) || missing(downwind_control_subset) ||
        missing(positive_subset)) {
      stop("When eda_type = 'num_obs_days_by_year', you must supply: data, day_column_name, year_column_name, upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset, positive_subset.")
    }
  }

  if (eda_type == "hist_day_group_sizes") {
    if (missing(data) || missing(day_column_name) ||
        missing(upwind_subset) || missing(downwind_subset) || missing(positive_subset)) {
      stop("When eda_type = 'hist_day_group_sizes', you must supply: data, day_column_name, upwind_subset, downwind_subset, positive_subset.")
    }
  }

  if (eda_type == "qq_rain") {
    if (missing(data) || missing(rain_col_name) || missing(use_raw) ||
        missing(upwind_subset) || missing(downwind_subset) ||
        missing(downwind_target_subset) || missing(downwind_control_subset) ||
        missing(positive_subset)) {
      stop("When eda_type = 'qq_rain', you must supply: data, rain_col_name, use_raw, upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset, positive_subset.")
    }
  }

  if (eda_type == "ts_by_type") {
    if (missing(data) || missing(rain_col_name) || missing(day_column_name) ||
        missing(year_column_name) || missing(use_raw) ||
        missing(upwind_subset) || missing(downwind_target_subset) || missing(downwind_control_subset) ||
        missing(positive_subset)) {
      stop("When eda_type = 'ts_by_type', you must supply: data, rain_col_name, day_column_name, year_column_name, use_raw, upwind_subset, downwind_target_subset, downwind_control_subset, positive_subset.")
    }
  }

  if (eda_type == "ts_by_gauge") {
    if (missing(data) || missing(rain_col_name) || missing(day_column_name) ||
        missing(year_column_name) || missing(use_raw) || missing(gauge_id_column_name)) {
      stop("When eda_type = 'ts_by_gauge', you must supply: data, rain_col_name, day_column_name, year_column_name, use_raw, gauge_id_column_name.")
    }
  }

  if (eda_type == "ts_by_gauge_interactive") {
    if (missing(data) || missing(rain_col_name) || missing(day_column_name) ||
        missing(year_column_name) || missing(use_raw) || missing(gauge_id_column_name)) {
      stop("When eda_type = 'ts_by_gauge_interactive', you must supply: data, rain_col_name, day_column_name, year_column_name, use_raw, gauge_id_column_name.")
    }
  }

  if (eda_type == "map_static") {
    if (missing(data) || missing(rain_col_name) || missing(day_column_name) ||
        missing(year_column_name) || missing(use_raw) ||
        missing(longlat_column_names) || missing(long_lim) || missing(lat_lim) ||
        missing(positive_subset)) {
      stop("When eda_type = 'map_static', you must supply: data, rain_col_name, day_column_name, year_column_name, use_raw, longlat_column_names, long_lim, lat_lim, positive_subset.")
    }
  }

  if (eda_type == "map_dynamic") {
    if (missing(data) || missing(rain_col_name) || missing(day_column_name) ||
        missing(year_column_name) || missing(use_raw) ||
        missing(longlat_column_names) || missing(long_lim) || missing(lat_lim)) {
      stop("When eda_type = 'map_dynamic', you must supply: data, rain_col_name, day_column_name, year_column_name, use_raw, longlat_column_names, long_lim, lat_lim.")
    }
  }


  if(eda_type == 'num_obs_days'){
    upwind_expr = rlang::enquo(upwind_subset)
    upwind = rlang::eval_tidy(upwind_expr, data = data)

    downwind_expr = rlang::enquo(downwind_subset)
    downwind = rlang::eval_tidy(downwind_expr, data = data)


    positive_expr = rlang::enquo(positive_subset)
    positive = rlang::eval_tidy(positive_expr, data = data)

    downwind_target_expr = rlang::enquo(downwind_target_subset)
    target = rlang::eval_tidy(downwind_target_expr, data = data)

    downwind_control_expr = rlang::enquo(downwind_control_subset)
    control = rlang::eval_tidy(downwind_control_expr, data = data)

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
    upwind_expr = rlang::enquo(upwind_subset)
    upwind = rlang::eval_tidy(upwind_expr, data = data)

    downwind_expr = rlang::enquo(downwind_subset)
    downwind = rlang::eval_tidy(downwind_expr, data = data)


    positive_expr = rlang::enquo(positive_subset)
    positive = rlang::eval_tidy(positive_expr, data = data)

    downwind_target_expr = rlang::enquo(downwind_target_subset)
    target = rlang::eval_tidy(downwind_target_expr, data = data)

    downwind_control_expr = rlang::enquo(downwind_control_subset)
    control = rlang::eval_tidy(downwind_control_expr, data = data)

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
    upwind_expr = rlang::enquo(upwind_subset)
    upwind = rlang::eval_tidy(upwind_expr, data = data)

    downwind_expr = rlang::enquo(downwind_subset)
    downwind = rlang::eval_tidy(downwind_expr, data = data)


    positive_expr = rlang::enquo(positive_subset)
    positive = rlang::eval_tidy(positive_expr, data = data)


    day_values = data[,day_column_name]

    upwind_days = day_values[upwind & positive]
    upwind_sizes = as.numeric(table(upwind_days))

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
      ggplot2::geom_histogram(binwidth = 1, fill = 'grey', color='black') +
      ggplot2::labs(title = upwind_title, x = "Group Size", y = "Frequency") +
      ggplot2::theme_bw()

    p2 <- ggplot2::ggplot(data.frame(size = downwind_sizes), ggplot2::aes(x = size)) +
      ggplot2::geom_histogram(binwidth = 1, fill = 'grey', color='black') +
      ggplot2::labs(title = downwind_title, x = "Group Size", y = "Frequency") +
      ggplot2::theme_bw()


    return(list(
      upwind_positive_hist = p1,
      downwind_positive_hist = p2
    ))
  }

  if(eda_type == 'qq_rain'){
    upwind_expr = rlang::enquo(upwind_subset)
    upwind = rlang::eval_tidy(upwind_expr, data = data)

    downwind_expr = rlang::enquo(downwind_subset)
    downwind = rlang::eval_tidy(downwind_expr, data = data)


    positive_expr = rlang::enquo(positive_subset)
    positive = rlang::eval_tidy(positive_expr, data = data)

    downwind_target_expr = rlang::enquo(downwind_target_subset)
    target = rlang::eval_tidy(downwind_target_expr, data = data)

    downwind_control_expr = rlang::enquo(downwind_control_subset)
    control = rlang::eval_tidy(downwind_control_expr, data = data)

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


    return(list(
      upwind_positive_qq = p1,
      downwind_positive_qq = p2
    ))
  }

  if(eda_type == 'ts_by_type'){

    upwind_expr = rlang::enquo(upwind_subset)
    upwind = rlang::eval_tidy(upwind_expr, data = data)



    positive_expr = rlang::enquo(positive_subset)
    positive = rlang::eval_tidy(positive_expr, data = data)

    downwind_target_expr = rlang::enquo(downwind_target_subset)
    target = rlang::eval_tidy(downwind_target_expr, data = data)

    downwind_control_expr = rlang::enquo(downwind_control_subset)
    control = rlang::eval_tidy(downwind_control_expr, data = data)

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
        ggplot2::theme_bw() +
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



    return(output_plot)


  }

  if(eda_type == 'ts_by_gauge_interactive'){
    if(!use_raw){
      data[,rain_col_name] = log(data[,rain_col_name])
      rain_label = paste0("log(", rain_col_name, ")")
    }else{
      rain_label = rain_col_name
    }




    #browser()
    if(!is.null(ts_focus_gauge)){
      data[,gauge_id_column_name] <- factor(data[,gauge_id_column_name])
      output_plot = ggplot2::ggplot(data, ggplot2::aes(x = .data[[day_column_name]],
                                                       y = .data[[rain_col_name]],
                                                       group = .data[[gauge_id_column_name]],
                                                       text = paste("Gauge:", .data[[gauge_id_column_name]]))) +
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
        ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]]), scales = "free") +
        ggplot2::labs(title = paste0('Time Series Plots of ', rain_label, ' of All Gauges (Hover over points to see Gauge ID)'),
                      y = rain_label,
                      color = gauge_id_column_name) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom")
    }else{
      output_plot = ggplot2::ggplot(data, ggplot2::aes(x = .data[[day_column_name]],
                                                       y = .data[[rain_col_name]],
                                                       group = .data[[gauge_id_column_name]],
                                                       text = paste("Gauge:", .data[[gauge_id_column_name]]))) +
        ggplot2::geom_line(alpha = 0.8, color = "grey80") +
        ggplot2::geom_point(color = "grey80") +
        ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]]), scales = "free") +
        ggplot2::labs(title = paste0('Time Series Plots of ', rain_label, ' of All Gauges (Hover over points to see Gauge ID)'),
                      y = rain_label) +
        ggplot2::theme_bw()
    }


    plotly_output = plotly::layout(plotly::ggplotly(output_plot, tooltip = "text"),
                                   legend = list(
                                     orientation = "h",
                                     x = 0.5,
                                     xanchor = "center",
                                     y = -0.2
                                   ))



    return(plotly_output)
  }

  if(eda_type == "map_static"){

    positive_expr = rlang::enquo(positive_subset)
    positive = rlang::eval_tidy(positive_expr, data = data)


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







    if(!is.null(ionizer_location_df)){
      ionizer_long = reshape(
        ionizer_location_df,
        varying = setdiff(colnames(ionizer_location_df), c(ionizer_id_column_name, ionizer_longlat_column_names)),
        v.names = "Deployed",
        timevar = year_column_name,
        times = setdiff(colnames(ionizer_location_df), c(ionizer_id_column_name, ionizer_longlat_column_names)),
        direction = "long"
      )
      ionizer_long <- ionizer_long[ionizer_long$Deployed == 1,
                                   c(ionizer_id_column_name, ionizer_longlat_column_names, year_column_name)]
    }



    if(!is.null(input_sf)){
      points_sf <- sf::st_as_sf(
        positive_df_avg,
        coords = c(longlat_column_names[1], longlat_column_names[2]),
        crs = 4326
      )

      if(!is.null(ionizer_location_df)){
        ionizer_long_sf = sf::st_as_sf(
          ionizer_long,
          coords = c(ionizer_longlat_column_names[1], ionizer_longlat_column_names[2]),
          crs = 4326
        )
      }




      if(!is.null(ionizer_location_df)){
        output_plot = ggplot2::ggplot() +
          ggplot2::geom_sf(data = input_sf, fill = "white", color = "black") +
          ggplot2::geom_sf(data = points_sf, ggplot2::aes(color = .data[[rain_col_name]])) +
          ggplot2::geom_sf(data = ionizer_long_sf, shape = 8) +
          ggplot2::geom_text(
            data = ionizer_long,
            ggplot2::aes(x = .data[[ionizer_longlat_column_names[1]]],
                         y = .data[[ionizer_longlat_column_names[2]]], label = .data[[ionizer_id_column_name]]),
            vjust = -0.8,
            size = 3,
            fontface = "bold"
          ) +
          ggplot2::coord_sf(xlim = long_lim, ylim = lat_lim, expand = FALSE)+
          ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]])) +
          ggplot2::labs(title = paste0('Spatial Plots of ', rain_label, ' Averaged Across\nAll Days Satisfying ',
                                       deparse(original_args$positive_subset)),
                        x = 'Longitude', y = 'Latitude', color = rain_label) +
          ggplot2::scale_color_viridis_c(name = rain_label) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = 'bottom')
      }else{
        output_plot = ggplot2::ggplot() +
          ggplot2::geom_sf(data = input_sf, fill = "white", color = "black") +
          ggplot2::geom_sf(data = points_sf, ggplot2::aes(color = .data[[rain_col_name]])) +
          ggplot2::coord_sf(xlim = long_lim, ylim = lat_lim, expand = FALSE) +
          ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]])) +
          ggplot2::labs(title = paste0('Spatial Plots of ', rain_label, ' Averaged Across\nAll Days Satisfying ',
                                       deparse(original_args$positive_subset)),
                        x = 'Longitude', y = 'Latitude', color = rain_label) +
          ggplot2::scale_color_viridis_c() +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = 'bottom')
      }

    }else{

      if(!is.null(ionizer_location_df)){

        output_plot = ggplot2::ggplot(positive_df_avg, ggplot2::aes(x = .data[[longlat_column_names[1]]],
                                                                    y = .data[[longlat_column_names[2]]])) +
          ggplot2::borders()+
          ggplot2::xlim(long_lim) + ggplot2::ylim(lat_lim) +
          ggplot2::geom_point(ggplot2::aes(color = .data[[rain_col_name]])) +
          ggplot2::geom_point(data = ionizer_long, ggplot2::aes(x = .data[[ionizer_longlat_column_names[1]]],
                                                                y = .data[[ionizer_longlat_column_names[2]]]),
                              shape = 8) +
          ggplot2::geom_text(
            data = ionizer_long,
            ggplot2::aes(x = .data[[ionizer_longlat_column_names[1]]],
                         y = .data[[ionizer_longlat_column_names[2]]], label = .data[[ionizer_id_column_name]]),
            vjust = -0.8,
            size = 3,
            fontface = "bold"
          ) +
          ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]]), scales = "free") +
          ggplot2::labs(title = paste0('Spatial Plots of ', rain_label, ' Averaged Across\nAll Days Satisfying ',
                                       deparse(original_args$positive_subset)),
                        x = 'Longitude', y = 'Latitude', color = rain_label) +
          ggplot2::scale_color_viridis_c() +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = 'bottom')
      }else{

        output_plot = ggplot2::ggplot(positive_df_avg, ggplot2::aes(x = .data[[longlat_column_names[1]]],
                                                                    y = .data[[longlat_column_names[2]]])) +
          ggplot2::borders()+
          ggplot2::xlim(long_lim) + ggplot2::ylim(lat_lim) +
          ggplot2::geom_point(ggplot2::aes(color = .data[[rain_col_name]])) +
          ggplot2::facet_wrap(ggplot2::vars(.data[[year_column_name]]), scales = "free") +
          ggplot2::labs(title = paste0('Spatial Plots of ', rain_label, ' Averaged Across\nAll Days Satisfying ',
                                       deparse(original_args$positive_subset)),
                        x = 'Longitude', y = 'Latitude', color = rain_label) +
          ggplot2::scale_color_viridis_c() +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = 'bottom')
      }

    }

    if(elev_contour){
      bbox_sf <- sf::st_as_sf(
        sf::st_sfc(
          sf::st_polygon(list(rbind(
            c(long_lim[1], lat_lim[1]),  # bottom-left
            c(long_lim[1], lat_lim[2]),  # top-left
            c(long_lim[2], lat_lim[2]),  # top-right
            c(long_lim[2], lat_lim[1]),  # bottom-right
            c(long_lim[1], lat_lim[1])   # close polygon
          ))),
          crs = 4326
        )
      )

      elev_raster <- elevatr::get_elev_raster(
        locations = bbox_sf,
        z = elev_resolution,          # adjust zoom/resolution
        clip = "bbox",  # rectangular crop
      )

      elev_points <- raster::rasterToPoints(elev_raster)
      elev_df <- data.frame(
        x = elev_points[, 1],
        y = elev_points[, 2],
        elev = elev_points[, 3]
      )

      elev_df_land <- elev_df[elev_df$elev >= 0, ]


      output_plot = output_plot + ggplot2::geom_contour(
        data = elev_df_land,
        ggplot2::aes(
          x = x,
          y = y,
          z = elev
        ),
        linewidth = 0.1,
        color = "black",
        linetype = "dashed"
      )
    }

    return(output_plot)


  }

  if(eda_type == 'map_dynamic'){
    data[,'alpha'] = 1
    data[oman$Rain.Gauge.Measurement == 0,'alpha'] = 0.2

    if(!use_raw){
      data[,rain_col_name] = log(data[,rain_col_name])
      rain_label = paste0("log(", rain_col_name, ")")
    }else{
      rain_label = rain_col_name
    }

    data_df = data[, c(longlat_column_names, day_column_name, year_column_name, rain_col_name, 'alpha')]

    if(!is.null(ionizer_location_df)){
      ionizer_long = reshape(
        ionizer_location_df,
        varying = setdiff(colnames(ionizer_location_df), c(ionizer_id_column_name, ionizer_longlat_column_names)),
        v.names = "Deployed",
        timevar = year_column_name,
        times = setdiff(colnames(ionizer_location_df), c(ionizer_id_column_name, ionizer_longlat_column_names)),
        direction = "long"
      )
      ionizer_long <- ionizer_long[ionizer_long$Deployed == 1,
                                   c(ionizer_id_column_name, ionizer_longlat_column_names, year_column_name)]
    }


    if(elev_contour){
      bbox_sf <- sf::st_as_sf(
        sf::st_sfc(
          sf::st_polygon(list(rbind(
            c(long_lim[1], lat_lim[1]),  # bottom-left
            c(long_lim[1], lat_lim[2]),  # top-left
            c(long_lim[2], lat_lim[2]),  # top-right
            c(long_lim[2], lat_lim[1]),  # bottom-right
            c(long_lim[1], lat_lim[1])   # close polygon
          ))),
          crs = 4326
        )
      )

      elev_raster <- elevatr::get_elev_raster(
        locations = bbox_sf,
        z = elev_resolution,          # adjust zoom/resolution
        clip = "bbox",  # rectangular crop
      )

      elev_points <- raster::rasterToPoints(elev_raster)
      elev_df <- data.frame(
        x = elev_points[, 1],
        y = elev_points[, 2],
        elev = elev_points[, 3]
      )

      elev_df_land <- elev_df[elev_df$elev >= 0, ]

    }

    if(!is.null(input_sf)){
      points_sf <- sf::st_as_sf(
        data_df,
        coords = c(longlat_column_names[1], longlat_column_names[2]),
        crs = 4326
      )


      if(!is.null(focus_year)){
        points_sf_year <- points_sf[ points_sf[[year_column_name]] %in% focus_year, ]
        if(!is.null(ionizer_location_df)){
          ionizer_long_day <- merge(unique(data_df[data_df[[year_column_name]] %in% focus_year, c(day_column_name, year_column_name)]), ionizer_long, by = year_column_name, all.x = TRUE)
        }
      }else{
        points_sf_year <- points_sf
        if(!is.null(ionizer_location_df)){
          ionizer_long_day <- merge(unique(data_df[, c(day_column_name, year_column_name)]), ionizer_long, by = year_column_name, all.x = TRUE)
        }
      }

      if(!is.null(ionizer_location_df)){
        ionizer_long_day_sf = sf::st_as_sf(
          ionizer_long_day,
          coords = c(ionizer_longlat_column_names[1], ionizer_longlat_column_names[2]),
          crs = 4326
        )


        output_plot = ggplot2::ggplot() +
          ggplot2::geom_sf(data = input_sf, fill = "white", color = "black", linewidth = 0.3) +
          ggplot2::geom_sf(
            data = points_sf_year,
            ggplot2::aes(color = .data[[rain_col_name]], alpha = .data[['alpha']]),
            size = 2
          )   +
          ggplot2::geom_sf(data = ionizer_long_day_sf, shape = 8) +
          ggplot2::geom_text(
            data = ionizer_long_day,
            ggplot2::aes(x = .data[[ionizer_longlat_column_names[1]]],
                         y = .data[[ionizer_longlat_column_names[2]]], label = .data[[ionizer_id_column_name]]),
            vjust = -0.8,
            size = 3,
            fontface = "bold"
          )  +
          ggplot2::coord_sf(xlim = long_lim, ylim = lat_lim, expand = FALSE) +
          ggplot2::labs(
            title = paste0(
              rain_label, " for ", day_column_name, " : {closest_state}"
            ),
            x = "Longitude",
            y = "Latitude",
            color = rain_label
          ) +
          ggplot2::scale_color_viridis_c(name = rain_label) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::guides(alpha = 'none') +
          gganimate::transition_states(
            states = .data[[day_column_name]],
            state_length = 1,
            transition_length = 1
          ) +
          gganimate::ease_aes("linear")

        if(elev_contour){
          output_plot = output_plot + ggplot2::geom_contour(
            data = elev_df_land,
            ggplot2::aes(
              x = x,
              y = y,
              z = elev
            ),
            linewidth = 0.1,
            color = "black",
            linetype = "dashed"
          )
        }

        output_animate = gganimate::animate(output_plot, nframes = length(unique(points_sf_year[[day_column_name]])), fps = fps, end_pause = 20)

      }else{
        output_plot = ggplot2::ggplot() +
          ggplot2::geom_sf(data = input_sf, fill = "white", color = "black", linewidth = 0.3) +
          ggplot2::geom_sf(
            data = points_sf_year,
            ggplot2::aes(color = .data[[rain_col_name]], alpha = .data[['alpha']]),
            size = 2
          )  +
          ggplot2::coord_sf(xlim = long_lim, ylim = lat_lim, expand = FALSE) +
          ggplot2::labs(
            title = paste0(
              rain_label, " for ", day_column_name, " : {closest_state}"
            ),
            x = "Longitude",
            y = "Latitude",
            color = rain_label
          ) +
          ggplot2::scale_color_viridis_c(name = rain_label) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::guides(alpha = 'none') +
          gganimate::transition_states(
            states = .data[[day_column_name]],
            state_length = 1,
            transition_length = 1
          ) +
          gganimate::ease_aes("linear")

        if(elev_contour){
          output_plot = output_plot + ggplot2::geom_contour(
            data = elev_df_land,
            ggplot2::aes(
              x = x,
              y = y,
              z = elev
            ),
            linewidth = 0.1,
            color = "black",
            linetype = "dashed"
          )
        }
        output_animate = gganimate::animate(output_plot, nframes = length(unique(points_sf_year[[day_column_name]])), fps = fps, end_pause = 20)
      }
    }else{
      if(!is.null(focus_year)){
        data_df_year <- data_df[ data_df[[year_column_name]] %in% focus_year, ]
        if(!is.null(ionizer_location_df)){
          ionizer_long_day <- merge(unique(data_df[data_df[[year_column_name]] %in% focus_year, c(day_column_name, year_column_name)]), ionizer_long, by = year_column_name, all.x = TRUE)
        }
      }else{
        data_df_year <- data_df
        if(!is.null(ionizer_location_df)){
          ionizer_long_day <- merge(unique(data_df[, c(day_column_name, year_column_name)]), ionizer_long, by = year_column_name, all.x = TRUE)
        }
      }

      if(!is.null(ionizer_location_df)){
        output_plot = ggplot2::ggplot(data_df_year, ggplot2::aes(x = .data[[longlat_column_names[1]]],
                                                                 y = .data[[longlat_column_names[2]]])) +
          ggplot2::xlim(long_lim) + ggplot2::ylim(lat_lim) +
          ggplot2::geom_point(ggplot2::aes(color = .data[[rain_col_name]],
                                           alpha = .data[['alpha']])) +
          ggplot2::geom_point(data = ionizer_long_day, ggplot2::aes(x = .data[[ionizer_longlat_column_names[1]]],
                                                                    y = .data[[ionizer_longlat_column_names[2]]]),
                              shape = 8) +
          ggplot2::geom_text(
            data = ionizer_long_day,
            ggplot2::aes(x = .data[[ionizer_longlat_column_names[1]]],
                         y = .data[[ionizer_longlat_column_names[2]]], label = .data[[ionizer_id_column_name]]),
            vjust = -0.8,
            size = 3,
            fontface = "bold"
          ) +
          ggplot2::labs(
            title = paste0(rain_label, " for ", day_column_name, " : {closest_state}"),
            x = "Longitude",
            y = "Latitude",
            color = rain_label) +
          ggplot2::scale_color_viridis_c(name = rain_label) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = 'bottom')+
          ggplot2::guides(alpha = 'none') +
          gganimate::transition_states(
            states = .data[[day_column_name]],
            state_length = 1,
            transition_length = 1
          ) +
          gganimate::ease_aes("linear") +
          ggplot2::borders()

        if(elev_contour){
          output_plot = output_plot + ggplot2::geom_contour(
            data = elev_df_land,
            ggplot2::aes(
              x = x,
              y = y,
              z = elev
            ),
            linewidth = 0.1,
            color = "black",
            linetype = "dashed"
          )
        }

        output_animate = gganimate::animate(output_plot, nframes = length(unique(data_df_year[[day_column_name]])), fps = fps, end_pause = 20)
      }else{
        output_plot = ggplot2::ggplot(data_df_year, ggplot2::aes(x = .data[[longlat_column_names[1]]],
                                                                 y = .data[[longlat_column_names[2]]])) +
          ggplot2::xlim(long_lim) + ggplot2::ylim(lat_lim) +
          ggplot2::geom_point(ggplot2::aes(color = .data[[rain_col_name]],
                                           alpha = .data[['alpha']]))  +
          ggplot2::labs(
            title = paste0(rain_label, " for ", day_column_name, " : {closest_state}"),
            x = "Longitude",
            y = "Latitude",
            color = rain_label) +
          ggplot2::scale_color_viridis_c(name = rain_label) +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position = 'bottom')+
          ggplot2::guides(alpha = 'none') +
          gganimate::transition_states(
            states = .data[[day_column_name]],
            state_length = 1,
            transition_length = 1
          ) +
          gganimate::ease_aes("linear") +
          ggplot2::borders()

        if(elev_contour){
          output_plot = output_plot + ggplot2::geom_contour(
            data = elev_df_land,
            ggplot2::aes(
              x = x,
              y = y,
              z = elev
            ),
            linewidth = 0.1,
            color = "black",
            linetype = "dashed"
          )
        }

        output_animate = gganimate::animate(output_plot, nframes = length(unique(data_df_year[[day_column_name]])), fps = fps, end_pause = 20)
      }


    }

    if(!is.null(animate_filename)){
      gganimate::anim_save(filename = animate_filename, animation = output_animate)
    }

    return(list(gganim_object = output_plot,
                gif_image = output_animate))
  }

}




