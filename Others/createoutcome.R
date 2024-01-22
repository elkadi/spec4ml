#create an outcome column of dichotomous values of 0/1 based on the presence of value V in a column no ColoumNo
createoutcome<-function(M,V,ColoumNo) {
  M$Outcome<-rep(0,length(M))
  M[grep(V,M[[,ColoumNo]])]$Outcome<-1 #change
  M
}
