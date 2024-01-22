PLSCV<-function(M,V,ncompMax=20,validation="repeatedcv"){
  f<-paste0(V,"~spc")
  ctrl <- caret::trainControl(method=validation,repeats = 3)
  mod <- caret::train(as.formula(f),data=M,method = "pls",tuneGrid=expand.grid(ncomp=1:ncompMax), trControl=ctrl,metric="Rsquared")
  R2<-mod$results$Rsquared[as.numeric(mod$bestTune)]
  NC<-as.numeric(mod$bestTune)
  c(R2,NC)
}
