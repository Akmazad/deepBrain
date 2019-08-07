dataDir="/short/yr31/aa7970/azData/DeepBrain/Data/"
input.pos.datFile="H3K27ac_rnaSeq.Pos.dat"
input.neg.datFile="H3K27ac_rnaSeq.Neg.dat"
output.pos.datFile="H3K27ac_rnaSeq.Pos.tfSpecific.dat"
output.neg.datFile="H3K27ac_rnaSeq.Neg.tfSpecific.dat"
output.combinedFile="H3K27ac_rnaSeq.Pos.Neg.tfSpecific.dat"
tfProfFile="UCSC_Encode_wgEncodeAwgTfbsUniform_metadata_690_TF_profiles.csv"
selectedTFs="tf.genes.expr_0.5_perc.csv"

library("data.table")
tfProf=read.csv(paste0(dataDir,tfProfFile),header=T)[,c(2,6)]
tfs=read.csv(paste0(dataDir,selectedTFs),header=T)
pos.dat=read.csv(paste0(dataDir,input.pos.datFile), sep="\t", header=T)
neg.dat=read.csv(paste0(dataDir,input.neg.datFile), sep="\t", header=T)

#only for one run
fwrite(rbind(pos.dat,neg.dat),paste0(dataDir,"H3K27ac_rnaSeq.Pos.Neg.dat"), row.names=F, quote=F, sep="\t")

processAtf <- function(aTF,dat,tfProf){
  print(aTF)
  aTF.Prof = paste0(as.vector(tfProf[as.character(tfProf$Factor)==aTF,2]),".narrowPeak.overlaps.bed")
  aTF.Prof.dat = as.data.frame(dat[,aTF.Prof])
  agg.val =  data.frame( x1 = apply(aTF.Prof.dat[1:ncol(aTF.Prof.dat)], 1, sum))
  colnames(agg.val)=aTF
  agg.val= as.data.frame(ifelse(agg.val>0,1,0))
}
new.pos.dat = pos.dat[,c(1:9)]
new.pos.dat = cbind(new.pos.dat, apply(as.matrix(tfs),1,FUN=function(x) processAtf(x,pos.dat,tfProf)))
write.csv(new.pos.dat,paste0(dataDir,output.pos.datFile), row.names=F, quote=F, sep="\t")

new.neg.dat = neg.dat[,c(1:9)]
new.neg.dat = cbind(new.neg.dat, apply(as.matrix(tfs),1,FUN=function(x) processAtf(x,neg.dat,tfProf)))
write.csv(new.neg.dat,paste0(dataDir,output.neg.datFile), row.names=F, quote=F, sep="\t")

#write (pos+neg) combined
write.csv(rbind(new.pos.dat,new.neg.dat),paste0(dataDir,output.combinedFile), row.names=F, quote=F, sep="\t")
