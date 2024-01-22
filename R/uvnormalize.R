#' unit vector normalization 
#'
#' @param x a hyperspec object.
#' @return a hyperspec object of unit vector normalized preprocessed spectra.
#' @examples uvnormalize(x)
uvnormalize<- function(x) {x / sqrt(sum(x^2))}
