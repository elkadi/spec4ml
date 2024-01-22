RegCV<-function(M,V,comp=20,model= "pls",validation="repeatedcv"){
  f<-paste0(V,"~spc")
  ctrl <- caret::trainControl(method=validation,repeats = 3)
  mod <- caret::train(as.formula(f),data=M,method = model,tuneLength=comp, trControl=ctrl)
  R2<-mod$results$Rsquared[as.numeric(mod$bestTune)]
  NC<-as.numeric(mod$bestTune)
  c(R2,NC)
}
