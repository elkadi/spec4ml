#create PS for the spectra
spectraPs<-function(M,DCN,skip,f,porp=0.8,...){
  Index<-1:nrow(M)
  Indexed<-cbind(Index,M[,1:DCN])
  Indexed<-Indexed[-(skip),]
  m.out <- matchit(formula(f),
                   data=Indexed,
                   distance='logit')
  skipsps<-rep(0,length(skip))
  ps<-c(skipsps,m.out$distance)
  names(ps)<-NULL
  ps
}
