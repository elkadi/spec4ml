#' create a list of all spectra after 53 diverse preprocessing
#'
#' @param rawspectra a hyperspec object.
#' @return A list of 53 preprocessed spectra.
#' @examples
#' SpectraPreProcess2(spectra)
SpectraPreProcess2<-function(rawspectra){
  uvnormalize<- function(x) {x / sqrt(sum(x^2))}
  ####SMs######
  Sm3<-apply(rawspectra,1,savgol,fl=3, forder = 2, dorder = 0);colnames(Sm3$spc)<-colnames(rawspectra$spc)
  Sm5<-apply(rawspectra,1,savgol,fl=5, forder = 2, dorder = 0);colnames(Sm5$spc)<-colnames(rawspectra$spc)
  Sm11<-apply(rawspectra,1,savgol,fl=11, forder = 2, dorder = 0);colnames(Sm11$spc)<-colnames(rawspectra$spc)
  ##Snv
  snvh<-function(spectra){
    spectra[[]]<-standardNormalVariate(spectra[[]])
    spectra}
  Snv<-snvh(rawspectra)
  Sm5Snv<-snvh(Sm5)
  Sm11Snv<-snvh(Sm11)
  Sm3Snv<-snvh(Sm3)
  ##EMSC
  emsch<-function(spectra){
    spectra[[]]<-EMSC(spectra[[]])$corrected
    spectra}
  Emsc<-emsch(rawspectra)
  Sm5Emsc<-emsch(Sm5)
  Sm11Emsc<-emsch(Sm11)
  Sm3Emsc<-emsch(Sm3)
  ##Baseline
  spacing=as.numeric(colnames(rawspectra$spc)[2])-as.numeric(colnames(rawspectra$spc)[1])
  bl<-(rawspectra-spc.rubberband(rawspectra,spline = FALSE))[,, min+spacing ~ max-spacing]
  Sm5bl<-(Sm5-spc.rubberband(Sm5,spline = FALSE))[,, min+spacing ~ max-spacing]
  Sm3bl<-(Sm3-spc.rubberband(Sm3,spline = FALSE))[,, min+spacing ~ max-spacing]
  Sm11bl<-(Sm11-spc.rubberband(Sm11,spline = FALSE))[,, min+spacing ~ max-spacing]
  Snvbl<-(Snv-spc.rubberband(Snv,spline = FALSE))[,, min+spacing ~ max-spacing]
  Sm5Snvbl<-(Sm5Snv-spc.rubberband(Sm5Snv,spline = FALSE))[,, min+spacing ~ max-spacing]
  Sm3Snvbl<-(Sm3Snv-spc.rubberband(Sm3Snv,spline = FALSE))[,, min+spacing ~ max-spacing]
  Emscbl<-(Emsc-spc.rubberband(Emsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  Sm5Emscbl<-(Sm5Emsc-spc.rubberband(Sm5Emsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  Sm11Emscbl<-(Sm11Emsc-spc.rubberband(Sm11Emsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  Sm3Emscbl<-(Sm3Emsc-spc.rubberband(Sm3Emsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  #Derivatives
  ##1stDev
  D1S9<-apply(rawspectra,1,savgol,fl=9, forder = 2, dorder = 1);colnames(D1S9$spc)<-colnames(rawspectra$spc)
  D1S3<-apply(rawspectra,1,savgol,fl=3, forder = 2, dorder = 1);colnames(D1S3$spc)<-colnames(rawspectra$spc)
  ##2nDev
  D2S3<-apply(rawspectra,1,savgol,fl=3, forder = 2, dorder = 2);colnames(D2S3$spc)<-colnames(rawspectra$spc)
  D2S11<-apply(rawspectra,1,savgol,fl=11, forder = 2, dorder = 2);colnames(D2S11$spc)<-colnames(rawspectra$spc)
  #Normalize
  Uvn<-apply(rawspectra,1,uvnormalize);colnames(Uvn$spc)<-colnames(rawspectra$spc)
  Sm5Uvn<-apply(Sm5,1,uvnormalize);colnames(Sm5Uvn$spc)<-colnames(rawspectra$spc)
  Sm3Uvn<-apply(Sm3,1,uvnormalize);colnames(Sm3Uvn$spc)<-colnames(rawspectra$spc)
  Sm11Uvn<-apply(Sm11,1,uvnormalize);colnames(Sm11Uvn$spc)<-colnames(rawspectra$spc)
  SnvUvn<-apply(Snv,1,uvnormalize);colnames(SnvUvn$spc)<-colnames(rawspectra$spc)
  Sm3SnvUvn<-apply(Sm3Snv,1,uvnormalize);colnames(Sm3SnvUvn$spc)<-colnames(rawspectra$spc)
  Sm5SnvUvn<-apply(Sm5Snv,1,uvnormalize);colnames(Sm5SnvUvn$spc)<-colnames(rawspectra$spc)
  Sm11SnvUvn<-apply(Sm11Snv,1,uvnormalize);colnames(Sm11SnvUvn$spc)<-colnames(rawspectra$spc)
  EmscUvn<-apply(Emsc,1,uvnormalize);colnames(EmscUvn$spc)<-colnames(rawspectra$spc)
  Sm3EmscUvn<-apply(Sm3Emsc,1,uvnormalize);colnames(Sm3EmscUvn$spc)<-colnames(rawspectra$spc)
  Sm5EmscUvn<-apply(Sm5Emsc,1,uvnormalize);colnames(Sm5EmscUvn$spc)<-colnames(rawspectra$spc)
  Sm11EmscUvn<-apply(Sm11Emsc,1,uvnormalize);colnames(Sm11EmscUvn$spc)<-colnames(rawspectra$spc)
  blUvn<-apply(bl,1,uvnormalize);colnames(blUvn$spc)<-colnames(bl$spc)
  Sm5blUvn<-apply(Sm5bl,1,uvnormalize);colnames(Sm5blUvn$spc)<-colnames(bl$spc)
  Sm3blUvn<-apply(Sm3bl,1,uvnormalize);colnames(Sm3blUvn$spc)<-colnames(bl$spc)
  Sm11blUvn<-apply(Sm11bl,1,uvnormalize);colnames(Sm11blUvn$spc)<-colnames(bl$spc)
  EmscblUvn<-apply(Emscbl,1,uvnormalize);colnames(EmscblUvn$spc)<-colnames(bl$spc)
  Sm5EmscblUvn<-apply(Sm5Emscbl,1,uvnormalize);colnames(Sm5EmscblUvn$spc)<-colnames(bl$spc)
  Sm3EmscblUvn<-apply(Sm3Emscbl,1,uvnormalize);colnames(Sm3EmscblUvn$spc)<-colnames(bl$spc)
  Sm11EmscblUvn<-apply(Sm11Emscbl,1,uvnormalize);colnames(Sm11EmscblUvn$spc)<-colnames(bl$spc)
  D1S3Uvn<-apply(D1S3,1,uvnormalize);colnames(D1S3Uvn$spc)<-colnames(rawspectra$spc)
  D1S9Uvn<-apply(D1S9,1,uvnormalize);colnames(D1S9Uvn$spc)<-colnames(rawspectra$spc)
  D2S11Uvn<-apply(D2S11,1,uvnormalize);colnames(D2S11Uvn$spc)<-colnames(rawspectra$spc)

  x<-c(rawspectra,Sm3,Sm5,Sm11,Snv,Sm5Snv,Sm11Snv,Emsc,Sm5Emsc,Sm11Emsc,bl,Sm5bl,Sm11bl,Snvbl,Sm5Snvbl,Emscbl,Sm5Emscbl,Sm11Emscbl,D1S9,D2S11,Uvn,Sm5Uvn,Sm11Uvn,SnvUvn,Sm5SnvUvn,Sm11SnvUvn,EmscUvn,Sm5EmscUvn,Sm11EmscUvn,blUvn,Sm5blUvn,EmscblUvn,Sm5EmscblUvn,Sm11EmscblUvn,D1S9Uvn,D2S11Uvn,Sm3Snv,Sm3Emsc,Sm3bl,Sm3Snvbl,Sm3Emscbl,D1S3,D2S3,Sm3Uvn,Sm3SnvUvn,Sm3EmscUvn,Sm3blUvn,Sm11blUvn,Sm3EmscblUvn,D1S3Uvn)
}
