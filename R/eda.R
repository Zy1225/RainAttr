
#TODO: Add the following types of EDA

# Table summarizing number of observations/unique groups of Upwind/Downwind/Target/Control vs. Positive/Zero

# Histogram of cluster sizes for Upwind+Positive / Downwind+Positive

# QQplots of raw/log rainfall for Upwind+Positive / Downwind+Positive

# 2-dimensional heatmap plot for raw/lograinfall (including those with zeros), with row = gauge, column = day

# Time series plot (separate for each year) with each line representing mean (with and without zeros) raw/log rainfall of Upwind/Target/Control

# Time series plot (separate for each year) with each line representing the OBSERVED raw/log rainfall of each gauge, with optional filtering of gauges

# A spatial plot of map (separate for each year) with dots colored based on the mean raw/log rainfall (averaged within each year, with and without zeros) of each gauge


eda = function(eda_type,
               data, rain_col_name, day_column_name, year_column_name, use_raw, filter_gauge = NULL,
               upwind_subset, downwind_subset, downwind_target_subset, downwind_control_subset, positive_subset){

}
