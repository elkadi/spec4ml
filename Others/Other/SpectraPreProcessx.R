#' create a list of all spectra after 74 diverse preprocessing
#'
#' @param rawspectra a hyperspec object.
#' @return A list of 74 preprocessed spectra.
#' @examples
#' SpectraPreProcessX(spectra)
SpectraPreProcessX<-function(rawspectra,SmLn=3,SmMn=5,SmHn=11){
  uvnormalize<- function(x) {x / sqrt(sum(x^2))}
  ####SMs######
  SmL<-apply(rawspectra,1,savgol,fl=SmLn, forder = 2, dorder = 0);colnames(SmL$spc)<-colnames(rawspectra$spc)
  SmM<-apply(rawspectra,1,savgol,fl=SmMn, forder = 2, dorder = 0);colnames(SmM$spc)<-colnames(rawspectra$spc)
  SmH<-apply(rawspectra,1,savgol,fl=SmHn, forder = 2, dorder = 0);colnames(SmH$spc)<-colnames(rawspectra$spc)
  ##Snv
  snvh<-function(spectra){
    spectra[[]]<-standardNormalVariate(spectra[[]])
    spectra}
  Snv<-snvh(rawspectra)
  SmMSnv<-snvh(SmM)
  SmHSnv<-snvh(SmH)
  SmLSnv<-snvh(SmL)
  ##EMSC
  emsch<-function(spectra){
    spectra[[]]<-EMSC(spectra[[]])$corrected
    spectra}
  Emsc<-emsch(rawspectra)
  SmMEmsc<-emsch(SmM)
  SmHEmsc<-emsch(SmH)
  SmLEmsc<-emsch(SmL)
  ##Baseline
  spacing=as.numeric(colnames(rawspectra$spc)[2])-as.numeric(colnames(rawspectra$spc)[1])
  bl<-(rawspectra-spc.rubberband(rawspectra,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmMbl<-(SmM-spc.rubberband(SmM,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmLbl<-(SmL-spc.rubberband(SmL,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmHbl<-(SmH-spc.rubberband(SmH,spline = FALSE))[,, min+spacing ~ max-spacing]
  Snvbl<-(Snv-spc.rubberband(Snv,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmMSnvbl<-(SmMSnv-spc.rubberband(SmMSnv,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmLSnvbl<-(SmLSnv-spc.rubberband(SmLSnv,spline = FALSE))[,, min+spacing ~ max-spacing]
  Emscbl<-(Emsc-spc.rubberband(Emsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmMEmscbl<-(SmMEmsc-spc.rubberband(SmMEmsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmHEmscbl<-(SmHEmsc-spc.rubberband(SmHEmsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmLEmscbl<-(SmLEmsc-spc.rubberband(SmLEmsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  #Derivatives
  ##1stDev
  D1S9<-apply(rawspectra,1,savgol,fl=9, forder = 2, dorder = 1);colnames(D1S9$spc)<-colnames(rawspectra$spc)
  D1S3<-apply(rawspectra,1,savgol,fl=3, forder = 2, dorder = 1);colnames(D1S3$spc)<-colnames(rawspectra$spc)
  SmMD1S9<-apply(SmM,1,savgol,fl=9, forder = 2, dorder = 1);colnames(SmMD1S9$spc)<-colnames(rawspectra$spc)
  SmMD1S3<-apply(SmM,1,savgol,fl=3, forder = 2, dorder = 1);colnames(SmMD1S3$spc)<-colnames(rawspectra$spc)
  SmLD1S9<-apply(SmL,1,savgol,fl=9, forder = 2, dorder = 1);colnames(SmLD1S9$spc)<-colnames(rawspectra$spc)
  SmLD1S3<-apply(SmL,1,savgol,fl=3, forder = 2, dorder = 1);colnames(SmLD1S3$spc)<-colnames(rawspectra$spc)
  SmHD1S9<-apply(SmH,1,savgol,fl=9, forder = 2, dorder = 1);colnames(SmHD1S9$spc)<-colnames(rawspectra$spc)
  SmHD1S3<-apply(SmH,1,savgol,fl=3, forder = 2, dorder = 1);colnames(SmHD1S3$spc)<-colnames(rawspectra$spc)
  ##2nDev
  D2S3<-apply(rawspectra,1,savgol,fl=3, forder = 2, dorder = 2);colnames(D2S3$spc)<-colnames(rawspectra$spc)
  D2S11<-apply(rawspectra,1,savgol,fl=11, forder = 2, dorder = 2);colnames(D2S11$spc)<-colnames(rawspectra$spc)
  SmLD2S3<-apply(SmL,1,savgol,fl=3, forder = 2, dorder = 2);colnames(SmLD2S3$spc)<-colnames(rawspectra$spc)
  SmLD2S11<-apply(SmL,1,savgol,fl=11, forder = 2, dorder = 2);colnames(SmLD2S11$spc)<-colnames(rawspectra$spc)
  SmMD2S3<-apply(SmM,1,savgol,fl=3, forder = 2, dorder = 2);colnames(SmMD2S3$spc)<-colnames(rawspectra$spc)
  SmMD2S11<-apply(SmM,1,savgol,fl=11, forder = 2, dorder = 2);colnames(SmMD2S11$spc)<-colnames(rawspectra$spc)
  SmHD2S3<-apply(SmH,1,savgol,fl=3, forder = 2, dorder = 2);colnames(SmHD2S3$spc)<-colnames(rawspectra$spc)
  SmHD2S11<-apply(SmH,1,savgol,fl=11, forder = 2, dorder = 2);colnames(SmHD2S11$spc)<-colnames(rawspectra$spc)
  #Normalize
  Uvn<-apply(rawspectra,1,uvnormalize);colnames(Uvn$spc)<-colnames(rawspectra$spc)
  SmMUvn<-apply(SmM,1,uvnormalize);colnames(SmMUvn$spc)<-colnames(rawspectra$spc)
  SmLUvn<-apply(SmL,1,uvnormalize);colnames(SmLUvn$spc)<-colnames(rawspectra$spc)
  SmHUvn<-apply(SmH,1,uvnormalize);colnames(SmHUvn$spc)<-colnames(rawspectra$spc)
  SnvUvn<-apply(Snv,1,uvnormalize);colnames(SnvUvn$spc)<-colnames(rawspectra$spc)
  SmLSnvUvn<-apply(SmLSnv,1,uvnormalize);colnames(SmLSnvUvn$spc)<-colnames(rawspectra$spc)
  SmMSnvUvn<-apply(SmMSnv,1,uvnormalize);colnames(SmMSnvUvn$spc)<-colnames(rawspectra$spc)
  SmHSnvUvn<-apply(SmHSnv,1,uvnormalize);colnames(SmHSnvUvn$spc)<-colnames(rawspectra$spc)
  EmscUvn<-apply(Emsc,1,uvnormalize);colnames(EmscUvn$spc)<-colnames(rawspectra$spc)
  SmLEmscUvn<-apply(SmLEmsc,1,uvnormalize);colnames(SmLEmscUvn$spc)<-colnames(rawspectra$spc)
  SmMEmscUvn<-apply(SmMEmsc,1,uvnormalize);colnames(SmMEmscUvn$spc)<-colnames(rawspectra$spc)
  SmHEmscUvn<-apply(SmHEmsc,1,uvnormalize);colnames(SmHEmscUvn$spc)<-colnames(rawspectra$spc)
  blUvn<-apply(bl,1,uvnormalize);colnames(blUvn$spc)<-colnames(bl$spc)
  SmMblUvn<-apply(SmMbl,1,uvnormalize);colnames(SmMblUvn$spc)<-colnames(bl$spc)
  SmLblUvn<-apply(SmLbl,1,uvnormalize);colnames(SmLblUvn$spc)<-colnames(bl$spc)
  SmHblUvn<-apply(SmHbl,1,uvnormalize);colnames(SmHblUvn$spc)<-colnames(bl$spc)
  EmscblUvn<-apply(Emscbl,1,uvnormalize);colnames(EmscblUvn$spc)<-colnames(bl$spc)
  SmMEmscblUvn<-apply(SmMEmscbl,1,uvnormalize);colnames(SmMEmscblUvn$spc)<-colnames(bl$spc)
  SmLEmscblUvn<-apply(SmLEmscbl,1,uvnormalize);colnames(SmLEmscblUvn$spc)<-colnames(bl$spc)
  SmHEmscblUvn<-apply(SmHEmscbl,1,uvnormalize);colnames(SmHEmscblUvn$spc)<-colnames(bl$spc)
  D1S3Uvn<-apply(D1S3,1,uvnormalize);colnames(D1S3Uvn$spc)<-colnames(rawspectra$spc)
  D1S9Uvn<-apply(D1S9,1,uvnormalize);colnames(D1S9Uvn$spc)<-colnames(rawspectra$spc)
  D2S11Uvn<-apply(D2S11,1,uvnormalize);colnames(D2S11Uvn$spc)<-colnames(rawspectra$spc)

  SmMD1S9Uvn<-apply(SmMD1S9,1,uvnormalize);colnames(SmMD1S9Uvn$spc)<-colnames(rawspectra$spc)
  SmMD1S3Uvn<-apply(SmMD1S3,1,uvnormalize);colnames(SmMD1S3Uvn$spc)<-colnames(rawspectra$spc)
  SmLD1S9Uvn<-apply(SmLD1S9,1,uvnormalize);colnames(SmLD1S9Uvn$spc)<-colnames(rawspectra$spc)
  SmLD1S3Uvn<-apply(SmLD1S3,1,uvnormalize);colnames(SmLD1S3Uvn$spc)<-colnames(rawspectra$spc)
  SmHD1S9Uvn<-apply(SmHD1S9,1,uvnormalize);colnames(SmHD1S9Uvn$spc)<-colnames(rawspectra$spc)
  SmHD1S3Uvn<-apply(SmHD1S3,1,uvnormalize);colnames(SmHD1S3Uvn$spc)<-colnames(rawspectra$spc)
  SmLD2S3Uvn<-apply(SmLD2S3,1,uvnormalize);colnames(SmLD2S3Uvn$spc)<-colnames(rawspectra$spc)
  SmLD2S11Uvn<-apply(SmLD2S11,1,uvnormalize);colnames(SmLD2S11Uvn$spc)<-colnames(rawspectra$spc)
  SmMD2S3Uvn<-apply(SmMD2S3,1,uvnormalize);colnames(SmMD2S3Uvn$spc)<-colnames(rawspectra$spc)
  SmMD2S11Uvn<-apply(SmMD2S11,1,uvnormalize);colnames(SmMD2S11Uvn$spc)<-colnames(rawspectra$spc)
  SmHD2S3Uvn<-apply(SmHD2S3,1,uvnormalize);colnames(SmHD2S3Uvn$spc)<-colnames(rawspectra$spc)
  SmHD2S11Uvn<-apply(SmHD2S11,1,uvnormalize);colnames(SmHD2S11Uvn$spc)<-colnames(rawspectra$spc)
  x<-c(rawspectra,SmL,SmM,SmH,Snv,SmMSnv,SmHSnv,Emsc,SmMEmsc,SmHEmsc,bl,SmMbl,SmHbl,Snvbl,SmMSnvbl,Emscbl,SmMEmscbl,SmHEmscbl,D1S9,D2S11,Uvn,SmMUvn,SmHUvn,SnvUvn,SmMSnvUvn,SmHSnvUvn,EmscUvn,SmMEmscUvn,SmHEmscUvn,blUvn,SmMblUvn,EmscblUvn,SmMEmscblUvn,SmHEmscblUvn,D1S9Uvn,D2S11Uvn,SmLSnv,SmLEmsc,SmLbl,SmLSnvbl,SmLEmscbl,D1S3,D2S3,SmLUvn,SmLSnvUvn,SmLEmscUvn,SmLblUvn,SmHblUvn,SmLEmscblUvn,D1S3Uvn,SmMD1S9,SmMD1S3,SmLD1S9,SmLD1S3,SmHD1S9,SmHD1S3,SmLD2S3,SmLD2S11,SmMD2S3,SmMD2S11,SmHD2S3,SmHD2S11,SmMD1S9Uvn,SmMD1S3Uvn,SmLD1S9Uvn,SmLD1S3Uvn,SmHD1S9Uvn,SmHD1S3Uvn,SmLD2S3Uvn,SmLD2S11Uvn,SmMD2S3Uvn,SmMD2S11Uvn,SmHD2S3Uvn,SmHD2S11Uvn)
}
