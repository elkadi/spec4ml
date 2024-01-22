#' create a list of all spectra after 66 diverse preprocessing
#'
#' @param rawspectra a hyperspec object of the spectra to be prepossessed.
#' @param SmLn Numeric. The filter length of the low level of smoothing to be attempted.
#' @param SmHn Numeric. The filter length of the high level of smoothing to be attempted.
#' @param ExtraDerivativeLSm Numeric. The filter length of the low level of smoothing to be attempted for both the 1st and 2nd derivatives in addition to the default values (3).
#' @param ExtraDerivativeHSm Numeric. The filter length of the high level of smoothing to be attempted for both the 1st and 2nd derivatives in addition to the default values (9,11).
#' @return A list of 66 preprocessed spectra and a list of their names.
#' @import hyperSpec
#' @import pracma
#' @examples
#' SpectraPreProcess(spectra,SmLn=3,SmMn=5,SmHn=11, ExtraDerivativeLSm=15, ExtraDerivativeHSm=17)
SpectraPreProcess<-function(rawspectra,SmLn=3,SmMn=5,SmHn=11, ExtraDerivativeLSm=15, ExtraDerivativeHSm=17){
  ####SMs######
  SmL<-apply(rawspectra,1,pracma::savgol,fl=SmLn, forder = 2, dorder = 0);colnames(SmL$spc)<-colnames(rawspectra$spc)
  SmM<-apply(rawspectra,1,pracma::savgol,fl=SmMn, forder = 2, dorder = 0);colnames(SmM$spc)<-colnames(rawspectra$spc)
  SmH<-apply(rawspectra,1,pracma::savgol,fl=SmHn, forder = 2, dorder = 0);colnames(SmH$spc)<-colnames(rawspectra$spc)
  ##Snv
  Snv<-biospec::snvh(rawspectra)
  SmMSnv<-biospec::snvh(SmM)
  SmHSnv<-biospec::snvh(SmH)
  SmLSnv<-biospec::snvh(SmL)
  ##Baseline
  spacing=as.numeric(colnames(rawspectra$spc)[2])-as.numeric(colnames(rawspectra$spc)[1])
  bl<-(rawspectra-hyperSpec::spc.rubberband(rawspectra,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmMbl<-(SmM-hyperSpec::spc.rubberband(SmM,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmLbl<-(SmL-hyperSpec::spc.rubberband(SmL,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmHbl<-(SmH-hyperSpec::spc.rubberband(SmH,spline = FALSE))[,, min+spacing ~ max-spacing]
  Snvbl<-(Snv-hyperSpec::spc.rubberband(Snv,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmMSnvbl<-(SmMSnv-hyperSpec::spc.rubberband(SmMSnv,spline = FALSE))[,, min+spacing ~ max-spacing]
  SmLSnvbl<-(SmLSnv-hyperSpec::spc.rubberband(SmLSnv,spline = FALSE))[,, min+spacing ~ max-spacing]
  #Derivatives
  ##1stDev
  D1S9<-apply(rawspectra,1,pracma::savgol,fl=9, forder = 2, dorder = 1);colnames(D1S9$spc)<-colnames(rawspectra$spc)
  D1S3<-apply(rawspectra,1,pracma::savgol,fl=3, forder = 2, dorder = 1);colnames(D1S3$spc)<-colnames(rawspectra$spc)
  SmMD1S9<-apply(SmM,1,pracma::savgol,fl=9, forder = 2, dorder = 1);colnames(SmMD1S9$spc)<-colnames(rawspectra$spc)
  SmMD1S3<-apply(SmM,1,pracma::savgol,fl=3, forder = 2, dorder = 1);colnames(SmMD1S3$spc)<-colnames(rawspectra$spc)
  SmLD1S9<-apply(SmL,1,pracma::savgol,fl=9, forder = 2, dorder = 1);colnames(SmLD1S9$spc)<-colnames(rawspectra$spc)
  SmLD1S3<-apply(SmL,1,pracma::savgol,fl=3, forder = 2, dorder = 1);colnames(SmLD1S3$spc)<-colnames(rawspectra$spc)
  SmHD1S9<-apply(SmH,1,pracma::savgol,fl=9, forder = 2, dorder = 1);colnames(SmHD1S9$spc)<-colnames(rawspectra$spc)
  SmHD1S3<-apply(SmH,1,pracma::savgol,fl=3, forder = 2, dorder = 1);colnames(SmHD1S3$spc)<-colnames(rawspectra$spc)

  D1SH<-apply(rawspectra,1,pracma::savgol,fl=ExtraDerivativeHSm, forder = 2, dorder = 1);colnames(D1SH$spc)<-colnames(rawspectra$spc)#added in the version 2 (24/5)
  D1SL<-apply(rawspectra,1,pracma::savgol,fl=ExtraDerivativeLSm, forder = 2, dorder = 1);colnames(D1SL$spc)<-colnames(rawspectra$spc)#added in the version 2 (24/5)
  ##2nDev
  D2S3<-apply(rawspectra,1,pracma::savgol,fl=3, forder = 2, dorder = 2);colnames(D2S3$spc)<-colnames(rawspectra$spc)
  D2S11<-apply(rawspectra,1,pracma::savgol,fl=11, forder = 2, dorder = 2);colnames(D2S11$spc)<-colnames(rawspectra$spc)
  SmLD2S3<-apply(SmL,1,pracma::savgol,fl=3, forder = 2, dorder = 2);colnames(SmLD2S3$spc)<-colnames(rawspectra$spc)
  SmLD2S11<-apply(SmL,1,pracma::savgol,fl=11, forder = 2, dorder = 2);colnames(SmLD2S11$spc)<-colnames(rawspectra$spc)
  SmMD2S3<-apply(SmM,1,pracma::savgol,fl=3, forder = 2, dorder = 2);colnames(SmMD2S3$spc)<-colnames(rawspectra$spc)
  SmMD2S11<-apply(SmM,1,pracma::savgol,fl=11, forder = 2, dorder = 2);colnames(SmMD2S11$spc)<-colnames(rawspectra$spc)
  SmHD2S3<-apply(SmH,1,pracma::savgol,fl=3, forder = 2, dorder = 2);colnames(SmHD2S3$spc)<-colnames(rawspectra$spc)
  SmHD2S11<-apply(SmH,1,pracma::savgol,fl=11, forder = 2, dorder = 2);colnames(SmHD2S11$spc)<-colnames(rawspectra$spc)

  D2SL<-apply(rawspectra,1,pracma::savgol,fl=ExtraDerivativeLSm, forder = 2, dorder = 2);colnames(D2SL$spc)<-colnames(rawspectra$spc)#added in the version 2 (24/5)
  D2SH<-apply(rawspectra,1,pracma::savgol,fl=ExtraDerivativeHSm, forder = 2, dorder = 2);colnames(D2SH$spc)<-colnames(rawspectra$spc)#added in the version 2 (24/5)
  #Normalize
  Uvn<-apply(rawspectra,1,biospec::uvnormalize);colnames(Uvn$spc)<-colnames(rawspectra$spc)
  SmMUvn<-apply(SmM,1,biospec::uvnormalize);colnames(SmMUvn$spc)<-colnames(rawspectra$spc)
  SmLUvn<-apply(SmL,1,biospec::uvnormalize);colnames(SmLUvn$spc)<-colnames(rawspectra$spc)
  SmHUvn<-apply(SmH,1,biospec::uvnormalize);colnames(SmHUvn$spc)<-colnames(rawspectra$spc)
  SnvUvn<-apply(Snv,1,biospec::uvnormalize);colnames(SnvUvn$spc)<-colnames(rawspectra$spc)
  SmLSnvUvn<-apply(SmLSnv,1,biospec::uvnormalize);colnames(SmLSnvUvn$spc)<-colnames(rawspectra$spc)
  SmMSnvUvn<-apply(SmMSnv,1,biospec::uvnormalize);colnames(SmMSnvUvn$spc)<-colnames(rawspectra$spc)
  SmHSnvUvn<-apply(SmHSnv,1,biospec::uvnormalize);colnames(SmHSnvUvn$spc)<-colnames(rawspectra$spc)
  blUvn<-apply(bl,1,biospec::uvnormalize);colnames(blUvn$spc)<-colnames(bl$spc)
  SmMblUvn<-apply(SmMbl,1,biospec::uvnormalize);colnames(SmMblUvn$spc)<-colnames(bl$spc)
  SmLblUvn<-apply(SmLbl,1,biospec::uvnormalize);colnames(SmLblUvn$spc)<-colnames(bl$spc)
  SmHblUvn<-apply(SmHbl,1,biospec::uvnormalize);colnames(SmHblUvn$spc)<-colnames(bl$spc)
  D1S3Uvn<-apply(D1S3,1,biospec::uvnormalize);colnames(D1S3Uvn$spc)<-colnames(rawspectra$spc)
  D1S9Uvn<-apply(D1S9,1,biospec::uvnormalize);colnames(D1S9Uvn$spc)<-colnames(rawspectra$spc)
  D2S11Uvn<-apply(D2S11,1,biospec::uvnormalize);colnames(D2S11Uvn$spc)<-colnames(rawspectra$spc)

  D1SHUvn<-apply(D1SH,1,biospec::uvnormalize);colnames(D1SHUvn$spc)<-colnames(rawspectra$spc)#added in the version 2 (24/5)
  D1SLUvn<-apply(D1SL,1,biospec::uvnormalize);colnames(D1SLUvn$spc)<-colnames(rawspectra$spc)#added in the version 2 (24/5)
  D2SHUvn<-apply(D2SH,1,biospec::uvnormalize);colnames(D2SHUvn$spc)<-colnames(rawspectra$spc)#added in the version 2 (24/5)
  D2SLUvn<-apply(D2SL,1,biospec::uvnormalize);colnames(D2SLUvn$spc)<-colnames(rawspectra$spc)#added in the version 2 (24/5)

  SmMD1S9Uvn<-apply(SmMD1S9,1,biospec::uvnormalize);colnames(SmMD1S9Uvn$spc)<-colnames(rawspectra$spc)
  SmMD1S3Uvn<-apply(SmMD1S3,1,biospec::uvnormalize);colnames(SmMD1S3Uvn$spc)<-colnames(rawspectra$spc)
  SmLD1S9Uvn<-apply(SmLD1S9,1,biospec::uvnormalize);colnames(SmLD1S9Uvn$spc)<-colnames(rawspectra$spc)
  SmLD1S3Uvn<-apply(SmLD1S3,1,biospec::uvnormalize);colnames(SmLD1S3Uvn$spc)<-colnames(rawspectra$spc)
  SmHD1S9Uvn<-apply(SmHD1S9,1,biospec::uvnormalize);colnames(SmHD1S9Uvn$spc)<-colnames(rawspectra$spc)
  SmHD1S3Uvn<-apply(SmHD1S3,1,biospec::uvnormalize);colnames(SmHD1S3Uvn$spc)<-colnames(rawspectra$spc)
  SmLD2S3Uvn<-apply(SmLD2S3,1,biospec::uvnormalize);colnames(SmLD2S3Uvn$spc)<-colnames(rawspectra$spc)
  SmLD2S11Uvn<-apply(SmLD2S11,1,biospec::uvnormalize);colnames(SmLD2S11Uvn$spc)<-colnames(rawspectra$spc)
  SmMD2S3Uvn<-apply(SmMD2S3,1,biospec::uvnormalize);colnames(SmMD2S3Uvn$spc)<-colnames(rawspectra$spc)
  SmMD2S11Uvn<-apply(SmMD2S11,1,biospec::uvnormalize);colnames(SmMD2S11Uvn$spc)<-colnames(rawspectra$spc)
  SmHD2S3Uvn<-apply(SmHD2S3,1,biospec::uvnormalize);colnames(SmHD2S3Uvn$spc)<-colnames(rawspectra$spc)
  SmHD2S11Uvn<-apply(SmHD2S11,1,biospec::uvnormalize);colnames(SmHD2S11Uvn$spc)<-colnames(rawspectra$spc)
  Ms<-c(rawspectra,SmL,SmM,SmH,Snv,SmMSnv,SmHSnv,bl,SmMbl,SmHbl,Snvbl,SmMSnvbl,D1S9,D2S11,Uvn,SmMUvn,SmHUvn,SnvUvn,SmMSnvUvn,SmHSnvUvn,blUvn,SmMblUvn,D1S9Uvn,D2S11Uvn,SmLSnv,SmLbl,SmLSnvbl,D1S3,D2S3,SmLUvn,SmLSnvUvn,SmLblUvn,SmHblUvn,D1S3Uvn,SmMD1S9,SmMD1S3,SmLD1S9,SmLD1S3,SmHD1S9,SmHD1S3,SmLD2S3,SmLD2S11,SmMD2S3,SmMD2S11,SmHD2S3,SmHD2S11,SmMD1S9Uvn,SmMD1S3Uvn,SmLD1S9Uvn,SmLD1S3Uvn,SmHD1S9Uvn,SmHD1S3Uvn,SmLD2S3Uvn,SmLD2S11Uvn,SmMD2S3Uvn,SmMD2S11Uvn,SmHD2S3Uvn,SmHD2S11Uvn,D1SH,D1SL,D2SH,D2SL,D1SHUvn,D1SLUvn,D2SHUvn,D2SLUvn)
  Mnames<-c("rawspectra","SmL","SmM","SmH","Snv","SmMSnv","SmHSnv","bl","SmMbl","SmHbl","Snvbl","SmMSnvbl","D1S9","D2S11","Uvn","SmMUvn","SmHUvn","SnvUvn","SmMSnvUvn","SmHSnvUvn","blUvn","SmMblUvn","D1S9Uvn","D2S11Uvn","SmLSnv","SmLbl","SmLSnvbl","D1S3","D2S3","SmLUvn","SmLSnvUvn","SmLblUvn","SmHblUvn","D1S3Uvn","SmMD1S9","SmMD1S3","SmLD1S9","SmLD1S3","SmHD1S9","SmHD1S3","SmLD2S3","SmLD2S11","SmMD2S3","SmMD2S11","SmHD2S3","SmHD2S11","SmMD1S9Uvn","SmMD1S3Uvn","SmLD1S9Uvn","SmLD1S3Uvn","SmHD1S9Uvn","SmHD1S3Uvn","SmLD2S3Uvn","SmLD2S11Uvn","SmMD2S3Uvn","SmMD2S11Uvn","SmHD2S3Uvn","SmHD2S11Uvn","D1SH","D1SL","D2SH","D2SL","D1SHUvn","D1SLUvn","D2SHUvn","D2SLUvn")
  list(Ms,Mnames)
  }
