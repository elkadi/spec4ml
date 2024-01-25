#' SpectraSmoothing
#'
#' create a list of all spectra after 74 diverse preprocessing
#'
#' @param rawspectra a hyperspec object.
#' @return A list of preprocessed spectra.
#' @import hyperSpec
#' @import pracma
#' @examples
#' SpectraSmoothing(rawspectra,filter_lengths=seq(3, 20, by = 2))

SpectraSmoothing <- function(rawspectra, filter_lengths=seq(3, 20, by = 2)) {
  SmL <- apply(rawspectra, 1, function(x) {
    apply(filter_lengths, 1, function(fl) {
      savgol(x, fl = fl, forder = 2, dorder = 0)
    })
  })

  colnames(SmL$spc) <- colnames(rawspectra$spc)
  return(SmL)
}

