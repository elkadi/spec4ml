\name{SpectraPreProcess2}
\alias{SpectraPreProcess2}
\title{Preprocess spectra with flexible SG smoothing windows}
\description{
This function preprocesses hyperspectral data using Savitzky–Golay (SG) smoothing,
derivative filtering, baseline correction, standard normal variate (SNV) normalisation
and unit vector normalisation.  Unlike the original \code{SpectraPreProcess}
function which uses fixed SG window lengths【591189474738661†L14-L18】,
\code{SpectraPreProcess2} applies every odd window length from 3 up to a
user-specified maximum for smoothing.  Helper functions and loops are used
to eliminate code duplication.
}
\usage{
SpectraPreProcess2(rawspectra, maxSmWin = 11, derivWins1 = c(3, 9),
  derivWins2 = c(3, 11), extraDerivSmWins = c(15, 17))
}
\arguments{
  \item{rawspectra}{A \code{hyperSpec} object containing raw spectra.}
  \item{maxSmWin}{Largest odd SG window length used for smoothing.  All odd lengths from 3 up to
    \code{maxSmWin} are applied.}
  \item{derivWins1}{First derivative SG window lengths.  Defaults to \code{c(3, 9)}.}
  \item{derivWins2}{Second derivative SG window lengths.  Defaults to \code{c(3, 11)}.}
  \item{extraDerivSmWins}{Additional derivative window lengths applied to both first and second derivatives【591189474738661†L44-L48】.}
}
\value{
A list containing two elements: a named list of processed \code{hyperSpec}
objects and a character vector of the corresponding names.
}
\examples{
## Not run:
# data(spectra)  # assume a hyperSpec object
# res <- SpectraPreProcess2(spectra, maxSmWin = 11)
# names(res[[1]])
## End(Not run)
}