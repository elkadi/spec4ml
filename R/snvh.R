#' snvh
#'
#' Preprocess spectra in a hyperspec format by SNV
#'
#' @param spectra a hyperspec object.
#' @return a hyperspec object of SNV preprocessed spectra.
#' @import prospectr
#' @export

snvh<-function(spectra){
  spectra[[]]<-prospectr::standardNormalVariate(spectra[[]])
  spectra}
