#' SpecImp
#'
#' read csv files separated with ";" in R
#'
#' @param x a csv file.
#' @param y the number of lines to be skipped.
#' @param delimiter A delimiter for specifying the boundary between columns. Defaults is ";"
#' @return a dataframe of the imported spectra.
#' @export
#' @importFrom utils read.csv
SpecImp <- function(x,y=8, delimiter=";") {
  w<-read.csv(x, sep = delimiter,header = FALSE, skip = y)
}
