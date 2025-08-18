\name{SpectraPreProcess4}
\alias{SpectraPreProcess4}
\title{Spectral preprocessing with MSC and EMSC corrections}
\description{
\code{SpectraPreProcess4} further extends \code{SpectraPreProcess3} by
incorporating Multiplicative Scatter Correction (MSC) and Extended
Multiplicative Signal Correction (EMSC).  MSC removes additive and
multiplicative scattering effects by aligning each spectrum to a
reference【189035542347124†screenshot】.  EMSC performs model-based background
correction and normalisation of spectra, handling scaling, polynomial baselines
and interferents【729034731224589†L19-L33】.  These corrections are applied to
raw and smoothed spectra in addition to the existing SG smoothing,
derivative, baseline and normalisation steps.
}
\usage{
SpectraPreProcess4(rawspectra, maxSmWin = 11,
  extraDerivSmWins = c(15, 17), ...)
}
\arguments{
  \item{rawspectra}{A \code{hyperSpec} object.}
  \item{maxSmWin}{Largest odd SG window length used for smoothing and derivative filters.}
  \item{extraDerivSmWins}{Optional vector of extra odd window lengths for derivative filters.}
  \item{...}{Additional arguments passed to \code{EMSC::EMSC} when computing EMSC corrections.}
}
\value{
A list with processed spectra and their names, including MSC- and
EMSC-corrected variants.
}
\examples{
## Not run:
# res <- SpectraPreProcess4(spectra, maxSmWin = 11)
## End(Not run)
}