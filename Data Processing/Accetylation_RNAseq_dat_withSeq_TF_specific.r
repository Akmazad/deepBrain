
dataDir="/short/yr31/aa7970/azData/DeepBrain/Data/"
tfProf=read.csv(paste0(dataDir,"UCSC_Encode_wgEncodeAwgTfbsUniform_metadata_690_TF_profiles.csv"),header=T)[,c(2,6)]
tfProf[,2]=paste0(tfProf[,2],".narrowPeak.overlaps.bed")
tfs=read.csv(paste0(dataDir,"tf.genes.expr_0.5_perc.csv"),header=F)
pos.dat=read.csv(paste0(dataDir,"test.pos.dat"), sep="\t", header=T)
processAtf <- function(aTF,dat,tfProf){
}

