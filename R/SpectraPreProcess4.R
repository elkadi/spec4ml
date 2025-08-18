#' SpectraPreProcessOptimized
#'
#' Preprocess a hyperspectral dataset using a flexible set of Savitzky–Golay (SG) smoothing windows
#' and a streamlined pipeline for common spectral corrections.  In contrast to the original
#' implementation provided in `spec4ml`, this version removes a large amount of code duplication
#' and exposes parameters that allow the user to process spectra with every odd‑sized SG window up
#' to a user specified maximum.  Additional derivatives and normalisation steps are still
#' supported via arguments but can easily be tailored for bespoke workflows.
#'  In addition to smoothing, derivatives, SNV, baseline and UV normalisation,
#'  this version now also applies Multiplicative Scatter Correction (MSC) using
#'  `prospectr::msc` and Extended Multiplicative Signal Correction (EMSC) using
#'  `EMSC::EMSC`.  MSC removes additive and multiplicative scattering effects by
#'  aligning each spectrum to a reference【189035542347124†screenshot】, while EMSC
#'  performs model-based correction for scaling, polynomial baselines and
#'  interferents【729034731224589†L19-L33】.
#'
#' @param rawspectra A `hyperSpec` object containing the raw spectra to preprocess.
#' @param maxSmWin Integer.  The largest Savitzky–Golay window length (must be odd) to be used for
#'   smoothing.  All odd window lengths from 3 up to this value will be applied.  Defaults to 11
#'   to match the original function.  If the provided value is even it will be reduced by one.
#' @param derivWins1 Integer vector.  First derivative SG window lengths.  Defaults to c(3, 9).
#'   Odd numbers larger than two should be supplied.  Additional values may be supplied via
#'   `extraDerivSmWins`.
#' @param derivWins2 Integer vector.  Second derivative SG window lengths.  Defaults to c(3, 11).
#'   Odd numbers larger than two should be supplied.  Additional values may be supplied via
#'   `extraDerivSmWins`.
#' @param extraDerivSmWins Integer vector.  Additional SG window lengths to apply to both first
#'   and second derivatives.  Defaults to c(15, 17) to mirror the optional derivative smoothing
#'   present in the original code【955124169944246†L44-L48】.
#' @return A list containing two elements: (1) a named list of processed `hyperSpec` objects and
#'   (2) a character vector of the corresponding names.  The output names reflect the processing
#'   steps applied, for example `Sm5` (smoothed with window 5) or `Sm7blUvn` (window 7 smoothed,
#'   baseline corrected and UV normalised).
#' @import hyperSpec
#' @import pracma
#' @import spec4ml
#' @import prospectr
#' @import EMSC
#' @export
SpectraPreProcessOptimized <- function(
  rawspectra,
  maxSmWin = 11,
  derivWins1 = NULL,
  derivWins2 = NULL,
  extraDerivSmWins = c(15, 17)
) {
  ## sanity checks ------------------------------------------------------------
  if (!inherits(rawspectra, "hyperSpec")) {
    stop("`rawspectra` must be a hyperSpec object")
  }
  # enforce odd maximum smoothing window
  if (maxSmWin %% 2 == 0) maxSmWin <- maxSmWin - 1
  # generate smoothing window lengths
  smoothWins <- seq(3, maxSmWin, by = 2)
  # derive derivative window lengths from smoothing range when not supplied
  if (is.null(derivWins1)) derivWins1 <- smoothWins
  if (is.null(derivWins2)) derivWins2 <- smoothWins
  # merge in any extra derivative window lengths specified by the user
  derivWins1 <- sort(unique(c(derivWins1, extraDerivSmWins)))
  derivWins2 <- sort(unique(c(derivWins2, extraDerivSmWins)))
  # ensure derivative windows are odd and >1
  derivWins1 <- derivWins1[derivWins1 %% 2 == 1 & derivWins1 > 1]
  derivWins2 <- derivWins2[derivWins2 %% 2 == 1 & derivWins2 > 1]

  ## helper functions ---------------------------------------------------------
  # Apply Savitzky–Golay filter with given window and derivative order
  sg_apply <- function(spectra, fl, dorder = 0) {
    res <- apply(spectra, 1, pracma::savgol, fl = fl, forder = 2, dorder = dorder)
    # copy wavelength column names from original spectra
    colnames(res$spc) <- colnames(spectra$spc)
    res
  }
  # Compute standard normal variate (SNV) normalisation
  snv_apply <- function(spectra) {
    res <- spec4ml::snvh(spectra)
    colnames(res$spc) <- colnames(spectra$spc)
    res
  }
  # Compute baseline corrected spectra using the rubberband method
  baseline_apply <- function(spectra) {
    spacing <- as.numeric(colnames(spectra$spc)[2]) - as.numeric(colnames(spectra$spc)[1])
    bl <- spectra - hyperSpec::spc.rubberband(spectra, spline = FALSE)
    # remove marginal points to avoid artefacts at the extremes
    bl[, , min + spacing ~ max - spacing]
  }
  # UV (unit vector) normalisation
  uvnorm_apply <- function(spectra) {
    res <- spec4ml::uvnormalize(spectra)
    colnames(res$spc) <- colnames(spectra$spc)
    res
  }

  # Multiplicative scatter correction (MSC) using prospectr::msc.  This method
  # aligns each spectrum to an ideal reference (default mean spectrum) and
  # corrects additive and multiplicative effects【189035542347124†screenshot】.
  msc_apply <- function(spectra, ref_spectrum = NULL) {
    # extract spectral matrix; prospectr::msc expects samples in rows
    X <- as.matrix(spectra$spc)
    res_matrix <- if (is.null(ref_spectrum)) {
      prospectr::msc(X)
    } else {
      prospectr::msc(X, ref_spectrum = ref_spectrum)
    }
    res <- spectra
    res$spc <- res_matrix
    colnames(res$spc) <- colnames(spectra$spc)
    res
  }

  # Extended multiplicative signal correction (EMSC).  EMSC is a model-based
  # correction that accounts for scaling, polynomial baselines and interferents【729034731224589†L19-L33】.
  # Here we use the default model provided by EMSC::EMSC.  Additional
  # parameters can be passed via ... if required.
  emsc_apply <- function(spectra, ...) {
    X <- as.matrix(spectra$spc)
    emsc_res <- EMSC::EMSC(X, ...)
    corrected <- emsc_res$corrected
    res <- spectra
    res$spc <- corrected
    colnames(res$spc) <- colnames(spectra$spc)
    res
  }

  ## smoothing ----------------------------------------------------------------
  # Generate smoothed spectra for each odd window length (computed earlier in `smoothWins`)
  smoothed <- lapply(smoothWins, function(fl) sg_apply(rawspectra, fl = fl, dorder = 0))
  names(smoothed) <- paste0("Sm", smoothWins)

  ## basic corrections on raw --------------------------------------------------
  # Standard normal variate and baseline correction on raw spectra
  Snv <- snv_apply(rawspectra)
  bl <- baseline_apply(rawspectra)
  SnvBl <- baseline_apply(Snv)
  # Unit vector normalisation
  Uvn <- uvnorm_apply(rawspectra)
  blUvn <- uvnorm_apply(bl)

  # Multiplicative scatter correction (MSC) and extended MSC (EMSC) on raw spectra
  MSC <- msc_apply(rawspectra)
  EMSC_basic <- emsc_apply(rawspectra)

  ## MSC and EMSC on smoothed spectra ----------------------------------------
  smoothedMSC <- lapply(smoothed, msc_apply)
  names(smoothedMSC) <- paste0(names(smoothed), "MSC")
  smoothedEMSC <- lapply(smoothed, emsc_apply)
  names(smoothedEMSC) <- paste0(names(smoothed), "EMSC")

  ## corrections on smoothed spectra ------------------------------------------
  # SNV on each smoothed spectrum
  smoothedSnv <- lapply(smoothed, snv_apply)
  names(smoothedSnv) <- paste0(names(smoothed), "Snv")
  # baseline on each smoothed spectrum
  smoothedBl <- lapply(smoothed, baseline_apply)
  names(smoothedBl) <- paste0(names(smoothed), "bl")
  # SNV then baseline (in that order) on each smoothed spectrum
  smoothedSnvBl <- mapply(function(snvSp, nm) {
    res <- baseline_apply(snvSp)
    res
  }, snvSp = smoothedSnv, nm = names(smoothedSnv), SIMPLIFY = FALSE)
  names(smoothedSnvBl) <- paste0(names(smoothed), "Snvbl")
  # UV normalisation on smoothed and baseline corrected spectra
  smoothedUvn <- lapply(smoothed, uvnorm_apply)
  names(smoothedUvn) <- paste0(names(smoothed), "Uvn")
  smoothedBlUvn <- lapply(smoothedBl, uvnorm_apply)
  names(smoothedBlUvn) <- paste0(names(smoothed), "blUvn")
  smoothedSnvUvn <- lapply(smoothedSnv, uvnorm_apply)
  names(smoothedSnvUvn) <- paste0(names(smoothed), "SnvUvn")
  smoothedSnvBlUvn <- lapply(smoothedSnvBl, uvnorm_apply)
  names(smoothedSnvBlUvn) <- paste0(names(smoothed), "SnvblUvn")

  ## derivatives ---------------------------------------------------------------
  # First derivative on raw and smoothed spectra
  d1Raw <- lapply(derivWins1, function(fl) sg_apply(rawspectra, fl = fl, dorder = 1))
  names(d1Raw) <- paste0("D1S", derivWins1)
  # Second derivative on raw and smoothed spectra
  d2Raw <- lapply(derivWins2, function(fl) sg_apply(rawspectra, fl = fl, dorder = 2))
  names(d2Raw) <- paste0("D2S", derivWins2)
  # First and second derivatives on each smoothed spectrum
  d1Smoothed <- list()
  d2Smoothed <- list()
  for (flSm in smoothWins) {
    sName <- paste0("Sm", flSm)
    spc <- smoothed[[sName]]
    # compute first derivative for each derivative window
    for (fl in derivWins1) {
      key <- paste0(sName, "D1S", fl)
      d1Smoothed[[key]] <- sg_apply(spc, fl = fl, dorder = 1)
    }
    # compute second derivative for each derivative window
    for (fl in derivWins2) {
      key <- paste0(sName, "D2S", fl)
      d2Smoothed[[key]] <- sg_apply(spc, fl = fl, dorder = 2)
    }
  }
  # Normalisation of derivatives
  d1RawUvn <- lapply(d1Raw, uvnorm_apply)
  names(d1RawUvn) <- paste0(names(d1Raw), "Uvn")
  d2RawUvn <- lapply(d2Raw, uvnorm_apply)
  names(d2RawUvn) <- paste0(names(d2Raw), "Uvn")
  d1SmoothedUvn <- lapply(d1Smoothed, uvnorm_apply)
  names(d1SmoothedUvn) <- paste0(names(d1Smoothed), "Uvn")
  d2SmoothedUvn <- lapply(d2Smoothed, uvnorm_apply)
  names(d2SmoothedUvn) <- paste0(names(d2Smoothed), "Uvn")

  ## collect all outputs -------------------------------------------------------
  # Start with raw spectra and simple corrections
  Ms <- list(
    rawspectra = rawspectra,
    Snv = Snv,
    bl = bl,
    Snvbl = SnvBl,
    Uvn = Uvn,
    blUvn = blUvn
    , MSC = MSC
    , EMSC = EMSC_basic
  )
  # append smoothed and their variants
  Ms <- c(Ms, smoothed, smoothedSnv, smoothedBl, smoothedSnvBl, smoothedUvn,
          smoothedBlUvn, smoothedSnvUvn, smoothedSnvBlUvn)
  # append MSC and EMSC for smoothed spectra
  Ms <- c(Ms, smoothedMSC, smoothedEMSC)
  # append derivative spectra
  Ms <- c(Ms, d1Raw, d2Raw, d1Smoothed, d2Smoothed,
          d1RawUvn, d2RawUvn, d1SmoothedUvn, d2SmoothedUvn)
  Mnames <- names(Ms)
  list(Ms, Mnames)
}