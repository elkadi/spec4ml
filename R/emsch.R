#' preprocess spectra in hyperspec formate by emsc 
#'
#' @param spectra a hyperspec object.
#' @return a hyperspec object of emsc preprocessed spectra.
#' @examples emsch(spectra)
emsch<-function(spectra){
  spectra[[]]<-EMSC::EMSC(spectra[[]])$corrected
  spectra}
