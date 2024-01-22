#' Import and merge csv files from Avasoft software in a single csv file suitable for analysis
#'
#' @param folder a folder with the Avasoft csv files.
#' @param spectrometer nir (default) or visible.
#' @param saturation the value that should be assigned in case of saturation (default=5).
#' @param file the merged file.
#' @return a single csv file with all spectra merged and suitable for analysis.
#' @import hyperSpec
#' @examples

avasoft_importt<-function(folder=getwd(),spectrometer="nir",saturation=5, file="MergedSpectra.csv"){
if(spectrometer=="nir"){spectrometerID="*_1503184U2"}else{spectrometerID="*_1503137U2"}
##importing Files all files in a folder ending with the spectrometer ID
files <- Sys.glob(paste0(folder,"/",spectrometerID,".TXT"))
files <- files [seq (1, length (files))]
ImportedFiles<- lapply (files, SpecImp)

for (i in 1:length(ImportedFiles)){
  Abs<-NULL
  for (wl in 1: nrow(ImportedFiles[[i]])){
  if((ImportedFiles[[i]]$V2[wl]-ImportedFiles[[i]]$V3[wl])<0){
    Abs[wl]<-saturation
  }
  else{
    Abs[wl]<-(-log10(abs((ImportedFiles[[i]]$V2[wl]-ImportedFiles[[i]]$V3[wl])/(ImportedFiles[[i]]$V4[wl]-ImportedFiles[[i]]$V3[wl]))))
  }
  }
  ImportedFiles[[i]]<-cbind(ImportedFiles[[i]],Abs)
}
for (i in 1:length(ImportedFiles)){
  ImportedFiles[[i]]<-ImportedFiles[[i]][-c(2,3,4)]
}
##Approbriate naming
names(ImportedFiles)<-list.files(folder, pattern = paste0(spectrometerID,".TXT"))
names(ImportedFiles)<-sub(".TXT", "", names(ImportedFiles))
names(ImportedFiles)<-gsub(spectrometerID, "", names(ImportedFiles))

WL<-ImportedFiles[[1]][1]
colnames(WL)<-"WL"
SampleNames<-names(ImportedFiles)
for (i in seq(1,length(files), by=1)){
  colnames(ImportedFiles[[i]])<-c("WL",SampleNames[i])
}
df<-WL
#For loop to merge the spectra
for (i in seq(1,length(files), by=1)){
  Sample<-names(ImportedFiles)[i]
  colnames(ImportedFiles[[i]][2])<-Sample
  df<-cbind(df,ImportedFiles[[i]][2])
}
tdf<-tnName(df)
write.csv(tdf, file,row.names = TRUE)
}
