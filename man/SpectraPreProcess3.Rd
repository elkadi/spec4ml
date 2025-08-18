\name{SpectraPreProcess3}
\alias{SpectraPreProcess3}
\title{Spectral preprocessing with matching derivative and smoothing windows}
\description{
\code{SpectraPreProcess3} extends \code{SpectraPreProcess2} by using the same set
of odd SG window lengths for both smoothing and derivative calculations.  First
and second derivatives are computed for every odd window up to the maximum
smoothing length, ensuring that derivatives and smoothed spectra are
comparable.  Extra derivative windows may still be supplied.
}
\usage{
SpectraPreProcess3(rawspectra, maxSmWin = 11, extraDerivSmWins = c(15, 17))
}
\arguments{
  \item{rawspectra}{A \code{hyperSpec} object.}
  \item{maxSmWin}{Largest odd SG window length used for smoothing and derivative filters.}
  \item{extraDerivSmWins}{Optional vector of odd window lengths added to the derivative filters.}
}
\value{
As in \code{SpectraPreProcess2}, returns a list containing the processed spectra
and their names.
}
\examples{
## Not run:
# res <- SpectraPreProcess3(spectra, maxSmWin = 13)
## End(Not run)
}