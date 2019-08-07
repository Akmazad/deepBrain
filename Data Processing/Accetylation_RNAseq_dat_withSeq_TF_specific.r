dataDir="/short/yr31/aa7970/azData/DeepBrain/Data/"
input.pos.datFile="test.pos.dat"
input.neg.datFile="test.neg.dat"
output.pos.datFile="pos.tfSpecific.dat"
output.neg.datFile="neg.tfSpecific.dat"
output.combinedFile="pos.neg.tfSpecific.dat"
tfProfFile="UCSC_Encode_wgEncodeAwgTfbsUniform_metadata_690_TF_profiles.csv"
selectedTFs="tf.genes.expr_0.5_perc.csv"  # total 128 TFs

library("data.table")
tfProf=fread(paste0(dataDir,tfProfFile),header=T)[,c(2,6)]
tfs=fread(paste0(dataDir,selectedTFs),header=T)
pos.dat=fread(paste0(dataDir,input.pos.datFile), sep="\t", header=T)
neg.dat=fread(paste0(dataDir,input.neg.datFile), sep="\t", header=T)

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
fwrite(new.pos.dat,paste0(dataDir,output.pos.datFile), row.names=F, quote=F, sep="\t")

new.neg.dat = neg.dat[,c(1:9)]
new.neg.dat = cbind(new.neg.dat, apply(as.matrix(tfs),1,FUN=function(x) processAtf(x,neg.dat,tfProf)))
fwrite(new.neg.dat,paste0(dataDir,output.neg.datFile), row.names=F, quote=F, sep="\t")

#write (pos+neg) combined
fwrite(rbind(new.pos.dat,new.neg.dat),paste0(dataDir,output.combinedFile), row.names=F, quote=F, sep="\t")

