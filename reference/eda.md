# Exploratory Data Analysis for Rainfall Enhancement Trial Data

Performs various exploratory data analyses (EDA) on rainfall enhancement
trial data, including summaries of observation types and positive
rainfall events, time series plots, Q-Q plots, and static or animated
spatial maps.

## Usage

``` r
eda(
  eda_type,
  data,
  rain_col_name,
  day_column_name,
  year_column_name,
  use_raw,
  gauge_id_column_name,
  ts_focus_gauge = NULL,
  longlat_column_names,
  long_lim,
  lat_lim,
  input_sf = NULL,
  ionizer_location_df = NULL,
  ionizer_id_column_name = NULL,
  ionizer_longlat_column_names = NULL,
  elev_contour = TRUE,
  elev_resolution = 2,
  wind_direction_column_name = NULL,
  wind_arrow_long_lat = NULL,
  wind_speed_column_name = NULL,
  wind_speed_scaling = NULL,
  focus_year = NULL,
  fps = 10,
  animate_filename = NULL,
  upwind_subset,
  downwind_subset,
  downwind_target_subset,
  downwind_control_subset,
  positive_subset
)
```

## Arguments

- eda_type:

  A character string specifying the type of EDA to perform. Must be one
  of `"num_obs_days"`, `"num_obs_days_by_year"`,
  `"hist_day_group_sizes"`, `"qq_rain"`, `"ts_by_type"`,
  `"ts_by_gauge"`, `"ts_by_gauge_interactive"`, `"map_static"`,
  `"map_dynamic"`.

- data:

  A data frame containing the rainfall enhancement trial data.

- rain_col_name:

  A character string that refers to the column name of the raw scale
  rainfall in `data`.

- day_column_name:

  A character string that refers to the column name of day identifiers
  in `data`.

- year_column_name:

  A character string that refers to the column name of year identifiers
  in `data`.

- use_raw:

  Logical. If TRUE, EDA is performed on raw scale rainfall. if FALSE,
  EDA is performed on log-transformed rainfall.

- gauge_id_column_name:

  A character string that refers to the column name of gauge identifiers
  in `data`.

- ts_focus_gauge:

  An optional vector of gauge identifiers to highlight in time series
  plots. If `ts_focus_gauge` is not supplied, no highlighting is done
  and all gauges are plotted with the same color.

- longlat_column_names:

  A character vector of length 2 specifying the column names of
  longitude and latitude in `data`, for plotting spatial plots.

- long_lim:

  A numeric vector of length 2 specifying longitude limits for spatial
  plots.

- lat_lim:

  A numeric vector of length 2 specifying latitude limits for spatial
  plots.

- input_sf:

  An optional
  [`sf`](https://r-spatial.github.io/sf/reference/st_as_sf.html) object
  for polygon layers (e.g., country borders) in spatial plots. For
  example,
  `rnaturalearth::ne_countries(scale = "large", country = "Oman", returnclass = "sf")`.
  If `input_sf` is not supplied, then a default map is drawn using
  [`ggplot2::borders()`](https://ggplot2.tidyverse.org/reference/annotation_borders.html).

- ionizer_location_df:

  An optional data frame containing the locations of ionizers, to be
  used in the plotting of static and animated maps. If
  `ionizer_location_df` is not supplied, no ionizers are plotted. The
  data frame must contain columns specified by `ionizer_id_column_name`
  and `ionizer_longlat_column_names`. It should also include additional
  columns whose names correspond to the years in
  `unique(data[[year_column_name]])`, with binary values indicating
  whether each ionizer has been deployed in each year.

- ionizer_id_column_name:

  An optional character string that refers to the column name of ionizer
  identifiers in `ionizer_location_df`. This must be supplied if
  `ionizer_location_df` is supplied.

- ionizer_longlat_column_names:

  An optional character vector of length 2 specifying the column names
  of longitude and latitude in `ionizer_location_df`, for plotting
  ionizers. This must be supplied if `ionizer_location_df` is supplied.

- elev_contour:

  An optional logical. If `TRUE`, elevation contour lines (obtained from
  Amazon Web Services Terrain Tiles) are added to the spatial plots. If
  `FALSE` (default), elevation contour lines are not included.

- elev_resolution:

  An optional integer between 1 and 14, specifying the resolution of the
  elevation data obtained from Amazon Web Services Terrain Tiles via
  [`get_elev_raster`](https://rdrr.io/pkg/elevatr/man/get_elev_raster.html).
  Defaults to 2. Higher values indicate higher resolution; see the `z`
  argument in
  [`get_elev_raster`](https://rdrr.io/pkg/elevatr/man/get_elev_raster.html)
  for more details.

- wind_direction_column_name:

  An optional character string that refers to the column name of the
  daily wind direction in `data`. The wind direction should be expressed
  in degree, e.g., 90 represents easterly wind. If this is supplied,
  arrows representing daily wind direction are plotted on the animated
  map. It is recommended to supply `long_lim` and `lat_lim` so that the
  map is approximately square, ensuring that the wind arrow lengths
  appear visually consistent across directions.

- wind_arrow_long_lat:

  An optional numeric vector of length 2 specifying the longitude and
  latitude for the starting point of the wind arrows. If not supplied,
  it is set to be
  `c( long_lim[2] - (long_lim[2] - long_lim[1]) / 5,lat_lim[2] - (lat_lim[2] - lat_lim[1]) / 5)`.
  This is only used when `wind_direction_column_name` is supplied.

- wind_speed_column_name:

  An optional character string that refers to the column name of the
  daily wind speed in `data`. If supplied, wind arrow lengths are equal
  to wind speed multiplied with `wind_speed_scaling`. Otherwise, arrow
  lengths are constant and equal to `wind_speed_scaling`. This is only
  used when `wind_direction_column_name` is supplied.

- wind_speed_scaling:

  An optional positive numeric value controlling arrow length. If
  `wind_speed_column_name` is supplied, wind arrow lengths are equal to
  wind speed multiplied with this value. Otherwise, all wind arrows have
  length equal to this value. Default is 0.1 when
  `wind_speed_column_name` is supplied, and 1 otherwise. This is only
  used when `wind_direction_column_name` is supplied.

- focus_year:

  An optional vector specifying years to filter for animated maps. If
  `focus_year` is not supplied, all years are included.

- fps:

  An optional numeric specifying frames per second for animated maps.
  Default is `10`.

- animate_filename:

  An optional character string specifying the file name used to save the
  animated map as a GIF. The file name should end with `".gif"`. If
  `animate_filename` is not supplied, no GIF is saved on disk.

- upwind_subset:

  A logical expression used to extract the relevant subset of
  observations from `data` to be used in the upwind (first stage) LMM
  fitting. For example, `Gauge.Day.Type == "Upwind"`.

- downwind_subset:

  A logical expression used to extract the relevant subset of
  observations from `data` to be used in the downwind (second stage) LMM
  fitting. For example, `Gauge.Day.Type %in% c("Target","Control")`.

- downwind_target_subset:

  A logical expression used to extract the relevant subset of downwind
  (second stage) observations from `data` that were exposed to treatment
  (operating ionizers). For example, `Gauge.Day.Type == "Target"`.

- downwind_control_subset:

  A logical expression used to extract the relevant subset of downwind
  (second stage) observations from `data` that were not exposed to
  treatment (operating ionizers). For example,
  `Gauge.Day.Type == "Control"`.

- positive_subset:

  A logical expression used to extract the relevant subset of
  observations from `data` with positive rainfall. For example,
  `Rain.Gauge.Measurement > 0`.

## Value

Depending on `eda_type`:

- Data summaries as lists (`"num_obs_days"`, `"num_obs_days_by_year"`)

- [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  objects (`"hist_day_group_sizes"`, `"qq_rain"`, `"ts_by_type"`,
  `"ts_by_gauge"`, `"map_static"`)

- Interactive [`plotly`](https://rdrr.io/pkg/plotly/man/plot_ly.html)
  object (`"ts_by_gauge_interactive"`)

- A list consisting of an animated
  [`gganim`](https://gganimate.com/reference/gganimate-package.html)
  object and its associated gif_image (`"map_dynamic"`)

## Details

Each `eda_type` has its own behavior and relevant arguments provided in
parentheses:

- `num_obs_days` (`data`, `day_column_name`, `upwind_subset`,
  `downwind_subset`, `downwind_target_subset`,
  `downwind_control_subset`, `positive_subset`):

  Computes contingency tables of the number of observations and number
  of unique days for each subset (`upwind_subset`, `downwind_subset`,
  `downwind_target_subset`, `downwind_control_subset`), with columns
  corresponding to positive vs zero rainfall events.

- `num_obs_days_by_year` (`data`, `day_column_name`, `year_column_name`,
  `upwind_subset`, `downwind_subset`, `downwind_target_subset`,
  `downwind_control_subset`, `positive_subset`):

  Same as `num_obs_days`, but computed separately for each year.

- `hist_day_group_sizes` (`data`, `day_column_name`, `upwind_subset`,
  `downwind_subset`, `positive_subset`):

  Plots histograms of group sizes for days in `upwind_subset` and
  `downwind_subset` with positive rainfall. The group size for a given
  day and type (`upwind_subset`, `downwind_subset`) is defined as the
  number of gauges satisfying that type and having positive rainfall on
  that day.

- `qq_rain` (`data`, `rain_col_name`, `use_raw`, `upwind_subset`,
  `downwind_subset`, `downwind_target_subset`,
  `downwind_control_subset`, `positive_subset`):

  Produces Normal Q-Q plots for rainfall values (raw or log-transformed)
  for all observations satisfying `upwind_subset & positive_subset`, as
  well as those observations satisfying
  `downwind_subset & positive_subset`.

- `ts_by_type` (`data`, `rain_col_name`, `day_column_name`,
  `year_column_name`, `use_raw`, `upwind_subset`,
  `downwind_target_subset`, `downwind_control_subset`,
  `positive_subset`):

  Plots daily average rainfall (raw or log-transformed) by subset
  (`upwind_subset`, `downwind_target_subset`,
  `downwind_control_subset`), facetted by year. Averaging is performed
  only over observations with positive rainfall within each day.

- `ts_by_gauge` (`data`, `rain_col_name`, `day_column_name`,
  `year_column_name`, `use_raw`, `gauge_id_column_name`,
  `ts_focus_gauge`):

  Plots daily rainfall (raw or log-transformed) time series for each
  gauge, optionally highlighting a subset of gauges, with faceting by
  year. When plotting log-transformed rainfall using
  [`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot.html),
  observations with zero rainfall are represented as a point at the
  bottommost of the plot.

- `ts_by_gauge_interactive` (`data`, `rain_col_name`, `day_column_name`,
  `year_column_name`, `use_raw`, `gauge_id_column_name`,
  `ts_focus_gauge`):

  Interactive version of `ts_by_gauge` using
  [`plotly`](https://rdrr.io/pkg/plotly/man/plot_ly.html). Points show
  gauge identifiers on hover. When plotting log-transformed rainfall
  (`use_raw = FALSE`), zero rainfall values result in `-Inf`. Unlike
  [`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot.html),
  [`plotly`](https://rdrr.io/pkg/plotly/man/plot_ly.html) will omit
  these points, breaking the lines. Therefore, users should be aware
  that lines may appear broken for days with zero rainfall when plotting
  log-transformed rainfall.

- `map_static` (`data`, `rain_col_name`, `day_column_name`,
  `year_column_name`, `use_raw`, `longlat_column_names`, `long_lim`,
  `lat_lim`, `input_sf`, `ionizer_location_df`,
  `ionizer_id_column_name`, `ionizer_longlat_column_names`,
  `elev_contour`, `elev_resolution`, `positive_subset`):

  Produces a static spatial map of annual average rainfall (raw or
  log-transformed), optionally overlaying an
  [`sf`](https://r-spatial.github.io/sf/reference/st_as_sf.html) polygon
  layer supplied via `input_sf` and adding elevation contour lines if
  `elev_contour = TRUE` as well as plotting ionizers supplied via
  `ionizer_location_df`, with faceting by year. Averaging is performed
  only over days with positive rainfall for each gauge. Requires the
  maps package to draw map borders when `input_sf` is not supplied. Also
  requires the elevatr and raster package to plot elevation contour
  lines (obtained from Amazon Web Services Terrain Tiles) when
  `elev_contour = TRUE`.

- `map_dynamic` (`data`, `rain_col_name`, `day_column_name`,
  `year_column_name`, `use_raw`, `longlat_column_names`, `long_lim`,
  `lat_lim`, `input_sf`, `ionizer_location_df`,
  `ionizer_id_column_name`, `ionizer_longlat_column_names`,
  `elev_contour`, `elev_resolution`, `wind_direction_column_name`,
  `wind_arrow_long_lat`, `wind_speed_column_name`, `wind_speed_scaling`,
  `focus_year`, `fps`, `animate_filename`):

  Produces an animated map showing rainfall (raw or log-transformed) for
  each day, optionally filtered by year using the argument `focus_year`
  and overlaying an
  [`sf`](https://r-spatial.github.io/sf/reference/st_as_sf.html) polygon
  layer supplied via `input_sf` as well as adding elevation contour
  lines if `elev_contour = TRUE`, plotting ionizers supplied via
  `ionizer_location_df`, and plotting wind arrows based on daily wind
  directions supplied via `wind_direction_column_name`. The animation
  frames are displayed in the order of `data[,day_column_name]`. Users
  should ensure that `data[,day_column_name]` contains values that can
  be meaningfully ordered (e.g., numeric or Date), rather than
  nominal/factor values, so the animation reflects the correct temporal
  progression. Requires the maps package to draw map borders when
  `input_sf` is not supplied and the gifski package for rendering
  animated map. Also requires the elevatr and raster package to plot
  elevation contour lines (obtained from Amazon Web Services Terrain
  Tiles) when `elev_contour = TRUE`. When `animate_filename` is
  supplied, the resulting gif_image of the animated map will be saved to
  the location specified by `animate_filename`.
