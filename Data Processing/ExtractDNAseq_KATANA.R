rm(list = ls(all.names = TRUE))
setwd('/srv/scratch/z3526914/DeepBrain/Data')
library(data.table)
library(dplyr)
# read binned DNA seq (human)
hg <- fread("hg19.binwise.fasta.200bp.bed", sep="\t", header=F)
colnames(hg) <- c("chr","start","end","dna.seq","id","strand")

# read non-zero bins
nonZerobins <- fread("HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.dat", sep="\t", header=T)
nonZerobins.seq <- cbind(nonZerobins,hg[which(hg$id %in% nonZerobins$id), "dna.seq"])
colnames(nonZerobins.seq) <- c(colnames(nonZerobins),"dna.seq")
fwrite(nonZerobins.seq, file="HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.dat", sep="\t", row.names=F, quote=F)