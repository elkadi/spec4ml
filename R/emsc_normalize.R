#' EMSC Normalization of NIR Spectra
#'
#' Applies Extended Multiplicative Signal Correction (EMSC) normalization
#' to Near-Infrared (NIR) spectral data using reference controls averaged
#' by a specified grouping variable (e.g., "Day").
#'
#' This function reads input spectral and control data, computes mean reference
#' spectra for each group, applies EMSC normalization to each group of spectra,
#' and optionally saves the normalized spectra to a CSV file.
#'
#' @param spectra_file Character string. Name of the CSV file containing spectra data.
#' @param controls_file Character string. Name of the CSV file containing control data.
#' @param output_file Character string. Name of the output CSV file to save normalized data.
#'   Defaults to `"Spectra_Normalized.csv"`.
#' @param group_col Character string. Name of the column used to group data for normalization
#'   (e.g., `"Day"`). Defaults to `"Day"`.
#' @param spectra_dir Character string. Directory path containing the input spectra and control files.
#'   Defaults to `"../InputSpectra"`.
#' @param output_dir Character string. Directory path to save the normalized output.
#'   Defaults to `"../Normalized_Spectra"`.
#' @param save_output Logical. If `TRUE`, the normalized spectra are saved to disk.
#'   If `FALSE`, only the resulting dataframe is returned. Defaults to `TRUE`.
#'
#' @return A data frame containing the EMSC-normalized spectra combined with metadata.
#'
#' @details
#' The function identifies the start of the spectral region by detecting
#' the first numeric column in the input data. Metadata columns are preserved
#' and reattached to the normalized spectra.
#'
#' For each group (as defined by `group_col`), the function:
#' 1. Computes the mean control spectrum (`NIRC_Averages`).
#' 2. Builds an EMSC model based on the mean spectrum.
#' 3. Normalizes each corresponding spectral subset using that model.
#'
#' @examples
#' \dontrun{
#' # Normalize and save to file
#' emsc_normalize(
#'   spectra_file = "NIR_spectra.csv",
#'   controls_file = "Control_spectra.csv"
#' )
#'
#' # Normalize and return data only
#' normalized_df <- emsc_normalize(
#'   spectra_file = "NIR_spectra.csv",
#'   controls_file = "Control_spectra.csv",
#'   save_output = FALSE
#' )
#' }
#'
#' @importFrom dplyr group_by summarize_all
#' @importFrom EMSC EMSC
#' @importFrom utils read.csv write.csv
#' @importFrom dplyr %>%
#' @importFrom rlang .data

#'
#' @export
emsc_normalize <- function(
    spectra_file,
    controls_file,
    output_file = "NIR_spectra.csv",
    group_col = "Day",
    spectra_dir = "../InputSpectra",
    output_dir = "../Normalized_Spectra",
    save_output = TRUE
) {

  # Read input data
  Spectra <- utils::read.csv(spectra_file, check.names = FALSE)
  Controls <- utils::read.csv(controls_file, check.names = FALSE)

  # Compute mean control spectra by group
  NIRC_Averages <- Controls %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarize_all(.funs = mean)

  # Identify first numeric (spectral) column
  mc_Spectra <- find_first_numeric_column(Spectra)
  mc_NIRC_Averages <- find_first_numeric_column(NIRC_Averages)

  # Split metadata and spectra
  Spectra_SPC <- Spectra[, mc_Spectra:ncol(Spectra)]
  Spectra_Meta <- Spectra[, 1:(mc_Spectra - 1)]

  # Initialize output data
  NIR_emsc_Normalized <- NULL

  # Loop through groups and apply EMSC normalization
  for (i in unique(NIRC_Averages[[group_col]])) {
    spec <- Spectra_SPC[Spectra_Meta[[group_col]] == i, ]
    ref <- NIRC_Averages[NIRC_Averages[[group_col]] == i, mc_NIRC_Averages:ncol(NIRC_Averages)]

    emscmodel <- EMSC::EMSC(ref)
    emscnormalizedx <- EMSC::EMSC(spec, model = emscmodel$model)

    emscnormalized <- cbind(
      Spectra_Meta[Spectra_Meta[[group_col]] == i, ],
      emscnormalizedx[[1]]
    )

    NIR_emsc_Normalized <- rbind(NIR_emsc_Normalized, emscnormalized)
  }

  # Save output (if requested)
  if (save_output) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message("📁 Created output directory: ", output_dir)
    }

    output_path <- file.path(output_dir, output_file)
    utils::write.csv(NIR_emsc_Normalized, output_path, row.names = FALSE)
    message("✅ Normalization complete. File saved to: ", output_path)
  } else {
    message("✅ Normalization complete. Data returned (not saved).")
  }

  return(NIR_emsc_Normalized)
}
