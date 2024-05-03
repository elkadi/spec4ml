#'Import_and_Merge_Spectra
#'
#' Import and merge csv files of spectral data in a single csv file suitable for analysis. Suitable for files with only 2 columns: 1 for the Wavelength/ Wavenumber and the other with the absorbance,transmittance or reflectance
#'
#' @param folder a folder with the txt/csv files with spectral data.
#' @param file the merged file.
#' @param file_extension The extension of the spectral data files. Default is ".CSV"
#' @param delimiter The delimiter used in the files with spectral data.
#' @param r_Skip number of lines to be skipped before start reading (skipping initial metadata for instance). Default is 0
#' @return a single csv file with all spectra merged and suitable for analysis.
#' @importFrom utils write.csv
#' @export

Import_and_Merge_Spectra<-function(folder=getwd(),file="MergedSpectra.csv",delimiter=",",file_extension=".CSV",r_Skip=0){

# Importing Files: all files in a folder ending with the specified file extension
files <- Sys.glob(file.path(folder, paste0("*", file_extension)))
ImportedFiles<- lapply (files, SpecImp,y=r_Skip,delimiter=delimiter)

# Check if files were imported
if (length(ImportedFiles) == 0) {
  stop("No files found with the specified extension.")
}

# Appropriate naming
# Sample names
sample_names <- tools::file_path_sans_ext(basename(files))
names(ImportedFiles) <- sample_names

###Wavelength names
WL<-ImportedFiles[[1]][1]
colnames(WL)<-"WL"

###Adding the sample names
for (i in seq_along(files)){
  colnames(ImportedFiles[[i]])<-c("WL",sample_names[i])
}
df<-WL
#For loop to merge the spectra
for (i in seq_along(files)){
  df<-cbind(df,ImportedFiles[[i]][2])
}
# Transpose and write to CSV
tdf<-tnName(df)
write.csv(tdf, file,row.names = TRUE)
}
