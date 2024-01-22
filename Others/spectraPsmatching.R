#perform PS matching for the spectra
spectraPsmatching<-function(j,M,DCN,skip,f,porp=0.8,...){
  set.seed(j)
  Index<-1:nrow(M)
  Indexed<-cbind(Index,M[,1:DCN])
  Indexed<-Indexed[-(skip),]
  m.out <- matchit(formula(f),
                   data=Indexed,
                   distance='logit')
  v1<-rownames(m.out$match.matrix)
  V1s<-sample(v1,porp*length(v1))
  V2s<-m.out$match.matrix[V1s,]
  Vxs<-as.numeric(c(1:12,V1s,V2s))
  Vxs<-na.omit(Vxs)
  c(j,Vxs)
}
