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
#' The width of the `tpcs` matrix MUST be greater than or equal to the length of tempaxis.
#'
#' @returns A matrix of averaged values, where columns are temperatures, and rows are individual tpc runs.
#' @export
#'
#' @examplesIf interactive()
#' nla_tpc_empirical()
nla_tpc_empirical <- function(tpcs, dailytemps, tempaxis = NULL) {
	# TODO: Add more input format checking.
	# Check input formats
	if(!inherits(tpcs, "matrix")) {
		warning("tpcs is not a matrix")
		tpcs_mat <- tryCatch({
			as.matrix(tpcs)
		}, error = function(e) {
			stop("tpcs is not a matrix and cannot be coerced into one!")
		})
	} else {
		tpcs_mat <- tpcs
	}

	if(!inherits(dailytemps, "matrix")) {
		warning("dailytemps is not a matrix")
		dt_mat <- tryCatch({
			as.matrix(dailytemps)
		}, error = function(e) {
			stop("dailytemps is not a matrix and cannot be coerced into one!")
		})
	} else {
		dt_mat <- dailytemps
	}

	# Also width of tpcs must be the same length as tempaxis
	if (is.null(tempaxis)) {
		if (!is.null(colnames(tpcs))) {
			tempaxis <- as.numeric(colnames(tpcs))
		} else {
			stop("tempaxis not provided, and no column names on tpcs to infer from!")
		}
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

#' @title Run rate summation on a TPC with a functional form
#'
#' @param tpc_func The tpc function to apply.
#' @param paramlist A named list of parameters (or vectors of parameters) for `tpc_func`. If `paramlist$labels` is provided, this will be used as the rownames of the resultant matrix.
#' @param dailytemps A matrix of temperature across a 24 hour period at various mean temperatures.
#' @param tempaxis_name Override the name of the temperature axis within the tpc func (if it is not "temp", the default in rTPC).
#' @param minvalue Minimum value a TPC can have (acts as a lower clamp).
#' @param maxvalue Maximum value a TPC can have (acts as a lower clamp).
#' @param default_value Value to replace NaN results from the TPC with, should these occur (this may be changed later to a possible function instead).
#'
#' @returns A matrix of averaged values, where columns are temperatures, and rows are individual tpc runs.
#'
#' @export
#'
#' @examples
#' test_param_list <- list(
#' "a" = c(0.01, 0.02),
#' "b" = 2,
#' "tmin" = 10,
#' "tmax" = 30,
#' "labels" = c("test", "test2")
#' )
#' test_dt <- dt_sin(dtr = 15, tbar = seq(20,45,0.1), sunrise = 6, resolution = 0.5, clamp = TRUE)
#' nla_tpc_analytic(
#' 	tpc_func = rTPC::briere2_1999,
#' 	paramlist = test_param_list,
#' 	dailytemps = test_dt,
#' 	minvalue = 0,
#' 	default_value = 0
#' 	)
#'

nla_tpc_analytic <- function(tpc_func, paramlist, dailytemps, tempaxis_name = "temp", minvalue = NULL, maxvalue = NULL, default_value = NULL) {
	# Check input formats
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
	colnames(dt_mat) <- NULL
	rownames(dt_mat) <- NULL

	# Extract labels and then remove these from the paramlist
	tpc_count <- max(sapply(paramlist, length))
	tpc_labels <- paramlist[["labels"]]
	paramlist[["labels"]] <- NULL
	# Check the formals of the tpc func to make sure only members are matches to tpc formals
	# Technically speaking, we could probably just ignore any that do not match, but we'll do this for now
	tpc_func_args <- formals(tpc_func)
	paramlist <- paramlist[names(paramlist)[names(paramlist) %in% names(tpc_func_args)]]
	# Now check for missing parameters
	# We could find all the formals with no default i.e inherits(formals(tpc_func)$param, "name")
	tpc_func_req_args <- names(tpc_func_args[sapply(tpc_func_args, inherits, what = "name")])
	# then check that they are in paramlist
	missing_req_args <- tpc_func_req_args[!(tpc_func_req_args %in% names(paramlist))]
	missing_req_args <- missing_req_args[missing_req_args != tempaxis_name]

	if (length(missing_req_args) > 0) {
		stop(paste("Missing required parameters: ", missing_req_args, collapse = ", "))
	}

	# unwrap dailytemps
	dt_vec <- as.vector(dailytemps)
	# Apply the function and then wrap back up to 3D as before
	# To do this we should probably convert the params to a data frame, which will fail if recycling is dodgy
	paramdf <- as.data.frame(paramlist)

	tpcout <- t(apply(paramdf, 1, \(x) {
		do.call(tpc_func, c(list(temp = dt_vec), as.list(x)))
	}))
	# Sort out value clamping to further restrict the output behaviour of TPCs
	if (!is.null(minvalue)) {
		tpcout <- pmax(tpcout, minvalue)
	}

	if (!is.null(maxvalue)) {
		tpcout <- pmin(tpcout, maxvalue)
	}

	if (!is.null(default_value)) {
		tpcout[which(is.na(tpcout))] <- default_value
	}

	# Reflow the matrix to make a 3D array
	dim(tpcout) <- c(nrow(tpcout), nrow(dt_mat), ncol(dt_mat))
	# Permute dimensions, so it goes timepoints, iterations, temps and then get the colMeans (average on first axis)
	final_array <- aperm(tpcout, c(2, 1, 3))
	final_out <- colMeans(final_array, dims = 1)

	# Reattach temp names
	colnames(final_out) <- dt_temp_names
	if (length(tpc_labels) == nrow(final_out)) {
		rownames(final_out) <- tpc_labels
	}

	return(final_out)
}
