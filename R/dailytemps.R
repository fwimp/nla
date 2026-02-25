#' Generate a daily temperature matrix based on a modified version of the Logan-Parton model.
#'
#' @param dtr Daily temperature range.
#' @param tbar A vector of mean temperatures of the day.
#' @param resolution Temporal resolution of the output (in terms of hours).
#' @param tempdecimalprecision Decimal precision of the resultant temperature.
#' @param clamp Whether to clamp the output values between the minimum and maximum values of `tbar`.
#'
#' @note If `length(tbar) == 1`, clamp is automatically set to `FALSE`.
#'
#' @section Modifications:
#' This is a heavily restricted and modified version of the Logan-Parton model. It is parameterised to a day length of 12 and a sunrise time of 6am.
#'
#' It makes sure that the mean and median of a given column equals the associated `tbar`
#'
#' @returns A matrix of `24/resolution` (rows) x `length(tbar)` (cols) containing temperature values at given times.
#'
#' @export
#'
#' @references Parton, W.J. and Logan, J.A., 1981. A model for diurnal variation in soil and air temperature. Agricultural meteorology, 23, pp.205-216.
#'
#' @examples
#' dt_modifiedlp(10, tbar = seq(20,45,0.1))
dt_modifiedlp <- function(dtr, tbar, resolution = 1, tempdecimalprecision = 1, clamp = TRUE) {
	# cli::cli_warn(
	# 	c("!" = "Warning: This has been modified from the original paper version",
	# 	  ">" = "FW (22/01/26): Added clamp parameter",
	# 	  ">" = "FW (23/01/26): Fixed rownames to be 0-24",
	# 	  ">" = "FW (23/01/26): Added resolution parameter"
	# 	),
	# 	.frequency = "regularly", .frequency_id = "LPcalc_modified")
	##################### Set parameters
	#      NOTE: Uses military time 0/24 = midnight; 1 = 1 AM; 14 = 2 PM, etc.
	#      NOTE: For the exp decay function to work, hours must continue to count up (past 23) and do not reset at midnight

	# Set parameters for day/night length, start times, and sin wave amplitude
	if (length(tbar) == 1) {
		clamp <- FALSE
	}
	day_length <- 12
	night_length <- 24 - day_length
	time_sunrise = 6
	time_sunset = time_sunrise + day_length
	amplitude <- dtr/2

	# Set other constants
	nocturnal_constant_tau = 4
	time_Tmax_afterNoon = 1.5
	correction_factor <- -0.0575824  # Coefficient for difference btw mean and median temp is specific for day_length = 12
	resdp <- decimalplaces(resolution)

	##################### Create hour & temperature data frames for doing math

	# Create vectors of hours based on sunrise/sunset times
	sin_hours <- seq(time_sunrise, (time_sunset - resolution)  , resolution) # hours that follow sin function (starts at sunrise)
	exp_hours <- seq(time_sunset,  (time_sunrise + (24 - resolution)), resolution) # hours that follow exp decay function (starts at sunset)

	# Create empty matrices for hours input
	sin_hour_input <- matrix(nrow = length(sin_hours), ncol = length(tbar))
	exp_hour_input <- matrix(nrow = length(exp_hours), ncol = length(tbar))

	# Fill them with the appropriate hours
	sin_hour_input[, 1:length(tbar)] <- sin_hours
	exp_hour_input[, 1:length(tbar)] <- exp_hours

	# Turn mean_temps into a matrix with one row, so we can rep() it below
	mean_temp_input <- matrix(nrow = 1, ncol = length(tbar))
	mean_temp_input[1, ] <- tbar

	# Create (filled) matrices for temp input
	sin_mean_temp_input <- mean_temp_input[rep(1, length(sin_hours)), ]
	exp_mean_temp_input <- mean_temp_input[rep(1, length(exp_hours)), ]

	# Create (filled) matrices for Tmedian, Tmin, Tmax, and Tsunset - need two versions/sizes for all except Tsunset
	sin_median_temp <- sin_mean_temp_input - correction_factor*dtr
	sin_Tmin <- sin_median_temp - amplitude
	sin_Tmax <- sin_median_temp + amplitude

	exp_median_temp <- exp_mean_temp_input - correction_factor*dtr
	exp_Tmin <- exp_median_temp - amplitude
	exp_Tmax <- exp_median_temp + amplitude
	exp_Tsunset = exp_Tmin + (exp_Tmax-exp_Tmin) * sin(pi*day_length / (day_length + 2*time_Tmax_afterNoon) )

	##################### Do the math & process the output

	# Calculate sin function with matrices
	sin_output <- sin_Tmin + (sin_Tmax-sin_Tmin) * sin(pi*(sin_hour_input - 12 + day_length/2) / (day_length + 2*time_Tmax_afterNoon))

	# Calculate exponential decay function with matrices
	exp_output <- (exp_Tmin - exp_Tsunset*exp(-night_length/nocturnal_constant_tau) + (exp_Tsunset-exp_Tmin)*exp(-(exp_hour_input-time_sunset)/nocturnal_constant_tau)) /
		(1 - exp(-night_length/nocturnal_constant_tau))

	# Combine sin and exponential function outputs into one matrix
	Logan_Parton_output <- rbind(sin_output, exp_output)

	# Round to nearest 0.1C
	Logan_Parton_output_rounded = round(Logan_Parton_output, tempdecimalprecision)

	if (clamp) {
		Logan_Parton_output_rounded[Logan_Parton_output_rounded < min(tbar)] <- min(tbar)
		Logan_Parton_output_rounded[Logan_Parton_output_rounded > max(tbar)] <- max(tbar)
	}

	# Add appropriate row and column names
	rownames(Logan_Parton_output_rounded) <- c(sin_hours, exp_hours) %% 24
	colnames(Logan_Parton_output_rounded) <- tbar

	# Return output
	return(Logan_Parton_output_rounded)

}

#' @title Generate a simple daily temperature sine curve
#'
#' @param dtr Daily temperature range.
#' @param tbar A vector of target mean temperatures.
#' @param sunrise The time that the sun rises.
#' @param resolution The resolution of the time series in hours.
#' @param clamp Whether to clamp the output between `min(tbar)` and `max(tbar)`, defaults to `TRUE`.
#'
#' @note If `length(tbar) == 1`, clamp is automatically set to `FALSE`.
#'
#' @returns A matrix of `24/resolution` (rows) x `length(tbar)` (cols) containing temperature values at given times.
#'
#' @export
#'
#' @examples
#' dt_sin(dtr = 15, tbar = seq(20,45,0.1), sunrise = 6, resolution = 0.5, clamp = TRUE)
dt_sin <- function(dtr, tbar, sunrise=6, resolution = 1, tempdecimalprecision = 1, clamp = TRUE) {
	if (length(tbar) == 1) {
		clamp <- FALSE
	}
	resdp <- decimalplaces(resolution)
	hourseq <- seq(0, 24 - resolution, resolution)
	hourseq_working <- hourseq - 6
	sinewave <- sin(((hourseq_working)/12) * pi)
	sindelta <- sinewave * (dtr/2)
	out<- sapply(tbar, \(x, sindelta) {
		x + sindelta
	},
	sindelta = sindelta)
	out <- round(out, tempdecimalprecision)

	if (clamp) {
		out[out < min(tbar)] <- min(tbar)
		out[out > max(tbar)] <- max(tbar)
	}

	rownames(out) <- round((hourseq + sunrise) %% 24, resdp)
	colnames(out) <- as.character(tbar)
	return(out)
}

#' @title Generate a simple daily temperature triangle curve
#'
#' @param dtr Daily temperature range.
#' @param tbar A vector of target mean temperatures.
#' @param sunrise The time that the sun rises.
#' @param resolution The resolution of the time series in hours.
#' @param clamp Whether to clamp the output between `min(tbar)` and `max(tbar)`, defaults to `TRUE`.
#'
#' @note If `length(tbar) == 1`, clamp is automatically set to `FALSE`.
#'
#' @returns A matrix of `24/resolution` (rows) x `length(tbar)` (cols) containing temperature values at given times.
#'
#' @export
#'
#' @examples
#' dt_tri(dtr = 15, tbar = seq(20,45,0.1), sunrise = 6, resolution = 0.5, clamp = TRUE)
dt_tri <- function(dtr, tbar, sunrise=6, resolution = 1, tempdecimalprecision = 1, clamp = TRUE) {
	if (length(tbar) == 1) {
		clamp <- FALSE
	}
	resdp <- decimalplaces(resolution)
	hourseq <- seq(0, 24 - resolution, resolution)
	hourseq_working <- hourseq

	triwave <- 2 * abs(2 * ((hourseq_working/24) - floor((hourseq_working/24) + 0.5))) - 1
	tridelta <- triwave * (dtr/2)

	out<- sapply(tbar, \(x, tridelta) {
		x + tridelta
	},
	tridelta = tridelta)
	# Round data to the same scale as the resolution
	out <- round(out, tempdecimalprecision)

	if (clamp) {
		out[out < min(tbar)] <- min(tbar)
		out[out > max(tbar)] <- max(tbar)
	}

	rownames(out) <- round((hourseq + sunrise) %% 24, resdp)
	colnames(out) <- as.character(tbar)
	return(out)
}
