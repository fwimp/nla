#' Perform nonlinear averaging of a temperature-dependent variable across a range of temperatures from empirically-derived data
#'
#' @param tpcs A matrix of values measured across a temperature axis (represented in `tempaxis`).
#' @param dailytemps A matrix of temperatures across a 24h period, equating to each entry in `tempaxis`.
#' @param tempaxis The basis temperature axis.
#'
#' @section Input dimensionality:
#'
#' This function relies on carefully and appropriately-formatted inputs.
#'
#' `tpcs` should be a 2D matrix where the rows correspond to individual tpcs or bootstrap iterations of a tpc, and columns correspond to the temperature axis in `tempaxis`.
#'
#' `dailytemps` should be a 2D matrix where the rows correspond to sampling times across a 24h window (so if your temperatures are hourly, this should have 24 rows), and columns correspond to mean temperatures. So if column 5 represents 10°C, `mean(dailytemps[,5]) == 10`
#'
#' @note
#' The precision of `tempaxis` MUST be greater than or equal to the values of `dailytemps`.
#'
#' The width of the `tpcs` matrix MUST be the same length as tempaxis.
#'
#' @returns A matrix of averaged values, where columns are temperatures, and rows are individual tpc runs.
#' @export
#'
#' @examplesIf interactive()
#' nla_tpc_empirical()
nla_tpc_empirical <- function(tpcs, dailytemps, tempaxis) {
	# Also width of tpcs must be the same length as tempaxis
	# Check input formats
	if(!inherits(tpcs, "matrix")) {
		warning("tpcs obj is not a matrix")
		tpcs_mat <- tryCatch({
			as.matrix(tpcs)
		}, error = function(e) {
			stop("tpcs obj is not a matrix and cannot be coerced into one")
		})
	} else {
		tpcs_mat <- tpcs
	}

	if(!inherits(dailytemps, "matrix")) {
		warning("dailytemps obj is not a matrix")
		dt_mat <- tryCatch({
			as.matrix(dailytemps)
		}, error = function(e) {
			stop("dailytemps obj is not a matrix and cannot be coerced into one")
		})
	} else {
		dt_mat <- dailytemps
	}

	# Save and remove names of dims
	dt_temp_names <- colnames(dt_mat)
	dt_time_names <- rownames(dt_mat)
	colnames(tpcs_mat) <- NULL
	rownames(tpcs_mat) <- NULL
	colnames(dt_mat) <- NULL
	rownames(dt_mat) <- NULL

	# Find decimal places to round dt_mat to, and do this
	decimal_places <- max(find_dps(tempaxis))
	dt_mat <- round(dt_mat, digits = decimal_places)
	# In order to match up the cell values with the tempaxis, we must convert both dt_mat and tempaxis to integers.
	# Note: at _extreme_ lengths of tempaxis or values in dt_mat, this method might break due to the integer limit. That should not be a concern though.
	dt_mat_ints <- as.integer(dt_mat * 10^decimal_places)
	tempaxis_ints <- as.integer(tempaxis * 10^decimal_places)

	# Actually match them up
	dt_mat_lookup <- match(dt_mat_ints, tempaxis_ints)

	# Lookup each cell as a column in tpcs to make a matrix of nrows(tpcs) x (nrow(dt_mat) * ncol(dt_mat))
	final_array <- tpcs_mat[,dt_mat_lookup]
	# Reflow the matrix to make a 3D array
	dim(final_array) <- c(nrow(tpcs_mat), nrow(dt_mat), ncol(dt_mat))
	# Permute dimensions, so it goes timepoints, iterations, temps and then get the colMeans (average on first axis)
	final_array <- aperm(final_array, c(2, 1, 3))
	final_out <- colMeans(final_array, dims = 1)

	# Reattach temp names
	colnames(final_out) <- dt_temp_names

	return(final_out)
}


