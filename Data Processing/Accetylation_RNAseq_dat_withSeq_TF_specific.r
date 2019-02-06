
dataDir="/short/yr31/aa7970/azData/DeepBrain/Data/"
tfProf=read.csv(paste0(dataDir,"UCSC_Encode_wgEncodeAwgTfbsUniform_metadata_690_TF_profiles.csv"),header=T)[,c(2,6)]
tfs=read.csv(paste0(dataDir,"tf.genes.expr_0.5_perc.csv"),header=T)
pos.dat=read.csv(paste0(dataDir,"test.pos.dat"), sep="\t", header=T)
neg.dat=read.csv(paste0(dataDir,"test.neg.dat"), sep="\t", header=T)

processAtf <- function(aTF,dat,tfProf){
  aTF.Prof = paste0(as.vector(tfProf[as.character(tfProf$Factor)==aTF,2]),".narrowPeak.overlaps.bed")
  aTF.Prof.dat = as.data.frame(dat[,aTF.Prof])
  agg.val =  data.frame( x1 = apply(aTF.Prof.dat[1:ncol(aTF.Prof.dat)], 1, sum))
  colnames(agg.val)=aTF
  agg.val= as.data.frame(ifelse(agg.val>0,1,0))
  #new.dat <<- cbind(new.dat,agg.val) ## assigning new values to the global variable 'new.pos.dat'
}
new.pos.dat = pos.dat[,c(1:9)]
new.pos.dat = cbind(new.pos.dat, apply(as.matrix(tfs),1,FUN=function(x) processAtf(x,pos.dat,tfProf)))
write.csv(new.dat,"pos.cl.csv")

new.neg.dat = neg.dat[,c(1:9)]
new.neg.dat = cbind(new.neg.dat, apply(as.matrix(tfs),1,FUN=function(x) processAtf(x,neg.dat,tfProf)))
write.csv(new.neg.dat,"neg.cl.csv")
