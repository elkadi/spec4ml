#' create a list of all spectra after 37 diverse preprocessing
#'
#' @param rawspectra a hyperspec object.
#' @return A list of 74 preprocessed spectra.
#' @examples
#' SpectraPreProcess3(spectra,SmLn=3,SmMn=5,SmHn=11)
SpectraPreProcess5<-function(rawspectra,SmLn=3){
  ####SMs######
  library(hyperSpec)
  SmL<-apply(rawspectra,1,pracma::savgol,fl=SmLn, forder = 2, dorder = 0);hyperSpec::colnames(SmL$spc)<-hyperSpec::colnames(rawspectra$spc)
  ##Snv
  Snv<-snvh(rawspectra)
  SmLSnv<-snvh(SmL)
  ##EMSC
  Emsc<-emsch(rawspectra)
  SmLEmsc<-emsch(SmL)
  ##Baseline
  spacing=as.numeric(hyperSpec::colnames(rawspectra$spc)[2])-as.numeric(hyperSpec::colnames(rawspectra$spc)[1])
  bl<-(rawspectra-hyperSpec::spc.rubberband(rawspectra,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmLbl<-(SmL-hyperSpec::spc.rubberband(SmL,spline = FALSE))[,, min+spacing ~ max-spacing]
  Snvbl<-(Snv-hyperSpec::spc.rubberband(Snv,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmLSnvbl<-(SmLSnv-hyperSpec::spc.rubberband(SmLSnv,spline = FALSE))[,, min+spacing ~ max-spacing]
  Emscbl<-(Emsc-hyperSpec::spc.rubberband(Emsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmLEmscbl<-(SmLEmsc-hyperSpec::spc.rubberband(SmLEmsc,spline = FALSE))[,, min+spacing ~ max-spacing]
  #Derivatives
  ##1stDev
  D1S9<-apply(rawspectra,1,pracma::savgol,fl=9, forder = 2, dorder = 1);hyperSpec::colnames(D1S9$spc)<-hyperSpec::colnames(rawspectra$spc)
  D1S3<-apply(rawspectra,1,pracma::savgol,fl=3, forder = 2, dorder = 1);hyperSpec::colnames(D1S3$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLD1S9<-apply(SmL,1,pracma::savgol,fl=9, forder = 2, dorder = 1);hyperSpec::colnames(SmLD1S9$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLD1S3<-apply(SmL,1,pracma::savgol,fl=3, forder = 2, dorder = 1);hyperSpec::colnames(SmLD1S3$spc)<-hyperSpec::colnames(rawspectra$spc)
  ##2nDev
  D2S3<-apply(rawspectra,1,pracma::savgol,fl=3, forder = 2, dorder = 2);hyperSpec::colnames(D2S3$spc)<-hyperSpec::colnames(rawspectra$spc)
  D2S11<-apply(rawspectra,1,pracma::savgol,fl=11, forder = 2, dorder = 2);hyperSpec::colnames(D2S11$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLD2S3<-apply(SmL,1,pracma::savgol,fl=3, forder = 2, dorder = 2);hyperSpec::colnames(SmLD2S3$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLD2S11<-apply(SmL,1,pracma::savgol,fl=11, forder = 2, dorder = 2);hyperSpec::colnames(SmLD2S11$spc)<-hyperSpec::colnames(rawspectra$spc)
  #Normalize
  Uvn<-apply(rawspectra,1,uvnormalize);hyperSpec::colnames(Uvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLUvn<-apply(SmL,1,uvnormalize);hyperSpec::colnames(SmLUvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  SnvUvn<-apply(Snv,1,uvnormalize);hyperSpec::colnames(SnvUvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLSnvUvn<-apply(SmLSnv,1,uvnormalize);hyperSpec::colnames(SmLSnvUvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  EmscUvn<-apply(Emsc,1,uvnormalize);hyperSpec::colnames(EmscUvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLEmscUvn<-apply(SmLEmsc,1,uvnormalize);hyperSpec::colnames(SmLEmscUvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  blUvn<-apply(bl,1,uvnormalize);hyperSpec::colnames(blUvn$spc)<-hyperSpec::colnames(bl$spc)
  SmLblUvn<-apply(SmLbl,1,uvnormalize);hyperSpec::colnames(SmLblUvn$spc)<-hyperSpec::colnames(bl$spc)
  EmscblUvn<-apply(Emscbl,1,uvnormalize);hyperSpec::colnames(EmscblUvn$spc)<-hyperSpec::colnames(bl$spc)
  SmLEmscblUvn<-apply(SmLEmscbl,1,uvnormalize);hyperSpec::colnames(SmLEmscblUvn$spc)<-hyperSpec::colnames(bl$spc)
  D1S3Uvn<-apply(D1S3,1,uvnormalize);hyperSpec::colnames(D1S3Uvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  D1S9Uvn<-apply(D1S9,1,uvnormalize);hyperSpec::colnames(D1S9Uvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  D2S11Uvn<-apply(D2S11,1,uvnormalize);hyperSpec::colnames(D2S11Uvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  
  SmLD1S9Uvn<-apply(SmLD1S9,1,uvnormalize);hyperSpec::colnames(SmLD1S9Uvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLD1S3Uvn<-apply(SmLD1S3,1,uvnormalize);hyperSpec::colnames(SmLD1S3Uvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLD2S3Uvn<-apply(SmLD2S3,1,uvnormalize);hyperSpec::colnames(SmLD2S3Uvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  SmLD2S11Uvn<-apply(SmLD2S11,1,uvnormalize);hyperSpec::colnames(SmLD2S11Uvn$spc)<-hyperSpec::colnames(rawspectra$spc)
  Ms<-c(rawspectra,SmL,Snv,Emsc,bl,Snvbl,Emscbl,D1S9,D2S11,Uvn,SnvUvn,EmscUvn,blUvn,EmscblUvn,D1S9Uvn,D2S11Uvn,SmLSnv,SmLEmsc,SmLbl,SmLSnvbl,SmLEmscbl,D1S3,D2S3,SmLUvn,SmLSnvUvn,SmLEmscUvn,SmLblUvn,SmLEmscblUvn,D1S3Uvn,SmLD1S9,SmLD1S3,SmLD2S3,SmLD2S11,SmLD1S9Uvn,SmLD1S3Uvn,SmLD2S3Uvn,SmLD2S11Uvn)
  Mnames<-c("rawspectra","SmL","Snv","Emsc","bl","Snvbl","Emscbl","D1S9","D2S11","Uvn","SnvUvn","EmscUvn","blUvn","EmscblUvn","D1S9Uvn","D2S11Uvn","SmLSnv","SmLEmsc","SmLbl","SmLSnvbl","SmLEmscbl","D1S3","D2S3","SmLUvn","SmLSnvUvn","SmLEmscUvn","SmLblUvn","SmLEmscblUvn","D1S3Uvn","SmLD1S9","SmLD1S3","SmLD2S3","SmLD2S11","SmLD1S9Uvn","SmLD1S3Uvn","SmLD2S3Uvn","SmLD2S11Uvn")
  list(Ms,Mnames)
}
