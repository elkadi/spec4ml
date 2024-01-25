#' RegCV
#'
#' Perform crossvalidation for regression models (This version need further revision)
#'
#' @param M a hyperspec object
#' @param V the traget variable
#' @param model the type of the model (default=pls).
#' @param validation the type of the CV (default=repeatedcv).
#' @return The R2 and and best tune parameters
#' @import hyperSpec
#' @examples RegCV(M,V,comp=20,model= "pls",validation="repeatedcv")

RegCV<-function(M,V,comp=20,model= "pls",validation="repeatedcv"){
  f<-paste0(V,"~spc")
  ctrl <- caret::trainControl(method=validation,repeats = 3)
  mod <- caret::train(as.formula(f),data=M,method = model,tuneLength=comp, trControl=ctrl)
  R2<-mod$results$Rsquared[as.numeric(mod$bestTune)]
  NC<-as.numeric(mod$bestTune)
  c(R2,NC)
}
