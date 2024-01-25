#' SpecImp
#'
#' read csv files separated with ";" in R
#'
#' @param x a csv file.
#' @param y the number of lines to be skipped.
#' @return a dataframe of the imported spectra.
#' @examples ImportedFiles<- lapply (files, SpecImp)
SpecImp <- function(x,y=8) {
  w<-read.csv(x, sep = ";",header = FALSE, skip = y)
}
