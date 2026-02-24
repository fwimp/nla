#' Find the decimal places in a vector of floats
#'
#' @param v A vector of floats.
#'
#' @returns A vector of the number of decimal places required to represent each float.
find_dps <- function(v) {
	sapply(v, decimalplaces)
}

#' Find the decimal places in a float
#'
#' @param x A float.
#'
#' @returns The number of decimal places required to represent the float.
decimalplaces <- function(x) {
	# Thanks to Gergely Daróczi on SO:
	# https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r
	if (abs(x - round(x)) > .Machine$double.eps^0.5) {
		# browser()
		xchar <- as.character(x)
		if (grepl("e-", xchar, fixed = TRUE)) {
			# if e- is in there use different approach
			as.numeric(strsplit(sub('0+$', '', xchar), "-", fixed = TRUE)[[1]][[2]])
		} else {
			nchar(strsplit(sub('0+$', '', xchar), ".", fixed = TRUE)[[1]][[2]])
		}
	} else {
		return(0)
	}
}
