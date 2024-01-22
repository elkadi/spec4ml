# This is a function that conduct PLSDA montocarlo double cross validation after PS matching
PLSDADCVApplyPs<-function(combinations,M,OutcomeIndex,USePs=TRUE){
  inTrain<-combinations[-1]
  training<-M[inTrain,]; testing<-M[-inTrain,]
  rownames(training$spc)<-rownames(training)
  rownames(testing$spc)<-rownames(testing)
  ctrl <- trainControl(classProbs = TRUE, summaryFunction = twoClassSummary)
  if (USePs==TRUE){X<-cbind(training$psscores,training$spc)
  testingX<-cbind(testing$psscores,testing$spc)
  } else
  {X<-training$spc
  testingX<-testing$spc}
  Y<-dummy_cols(as.matrix(training$Outcome))
  Y<-as.matrix(Y[-1])
  colnames(Y)<-c(0,1) #change based on the factors i.e 0,1 in case of using 0 instead of -1
  mod <- plsda(X,Y,trControl = ctrl, tuneGrid= expand.grid(ncomp=1:20))
  result<-predict(mod, newdata = testingX)
  xtab<-table(result, as.factor(testing$Outcome))
  seed<-combinations[1]
  Ac<-confusionMatrix(xtab,positive = "1")$overall[1]
  SPC<-mod$ncomp
  FP<-confusionMatrix(xtab,positive = "1")$table[2]
  FN<-confusionMatrix(xtab,positive = "1")$table[3]
  npos<-table(testing$Outcome)[[2]]
  nneg<-table(testing$Outcome)[[1]]
  SN<-as.numeric(confusionMatrix(xtab,positive = "1")$byClass[1])
  SP<-as.numeric(confusionMatrix(xtab,positive = "1")$byClass[2])
  cbind(seed,Ac,FP,FN,npos,nneg,SN,SP,SPC)
}
