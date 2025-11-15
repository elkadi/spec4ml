#' Identify the First Numeric (Spectral) Column
#'
#' Determines the index of the first numeric column in a data frame.
#' This function is useful for separating metadata columns from spectral
#' data in NIR or similar datasets, where spectral columns begin at the
#' first purely numeric column.
#'
#' @param data A data frame containing both metadata and spectral data.
#'
#' @return An integer giving the index of the first numeric column in the data frame.
#'
#' @details
#' The function checks each column name to determine whether it can be
#' converted to a numeric value. The first column name that successfully
#' converts to numeric (i.e., not `NA`) is considered the start of the spectral
#' region.
#'
#' For example, if your spectral columns are labeled by wavenumber (e.g., `1100`,
#' `1102`, `1104`, ...), this function returns the index of the first of those.
#'
#' @examples
#' \dontrun{
#' # Example dataframe
#' df <- data.frame(
#'   SampleID = c("A", "B", "C"),
#'   Day = c(1, 1, 2),
#'   `1100` = c(0.1, 0.2, 0.3),
#'   `1102` = c(0.2, 0.3, 0.4)
#' )
#'
#' # Find where spectra start
#' mc <- find_first_numeric_column(df)
#' print(mc)  # Should return 3
#' }
#'
#' @export
find_first_numeric_column <- function(data) {
  # Validate input
  if (!is.data.frame(data)) {
    stop("Input must be a data frame.")
  }
  
  # Identify first column name that can be converted to numeric
  idx <- which(sapply(names(data), function(x) !is.na(suppressWarnings(as.numeric(x)))))[1]
  
  if (is.na(idx)) {
    stop("No numeric column names found in the provided data frame.")
  }
  
  return(idx)
}
