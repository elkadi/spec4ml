# This is a function that produce statistics about the models of the PLSDADCV
modelperfromancesummary<-function(RM=RMpar){
  RM<-do.call(rbind.data.frame, RMpar)
  MPC<-as.numeric(names(sort(-table(RM$SPC)))[1])
  RMMPC<-RM[(RM$SPC==MPC),]
  ACCI<-conf.interval(RM$Ac)
  Total.Accuracy.ttest<-my.t.test(RM$Ac,mu = 0.5)
  cbind(V,MPC,round(mean(RM$Ac),5),paste(round(ACCI[[1]],3),"-",round(ACCI[[2]],3)), paste(round(min(RM$Ac),3),"-",round(max(RM$Ac),3)),mean(RM$AcPvalue),Total.Accuracy.ttest$p.value,mean(RM$Kp),as.numeric(names(sort(-table(RM$FP)))[1]),as.numeric(names(sort(-table(RM$FN)))[1]),round(mean(RM$SN),3),round(mean(RM$SP),3), RMMPC[order(-RMMPC$Ac),][1,1])
}
