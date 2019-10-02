rm(list = ls(all.names = TRUE))
setwd('/srv/scratch/z3526914/DeepBrain/Data')
library(data.table)
library(dplyr)

# ### Option list
args = commandArgs(trailingOnly=FALSE)

### Make Cluster nodes for parallelizing
dataDir = args[3]
flankingLength = args[4]
inputFile = args[5]
outputFile = args[5]

library("BSgenome.Hsapiens.UCSC.hg19")
hg <- BSgenome.Hsapiens.UCSC.hg19

# read non-zero bins
nonZerobins <- fread(inputFile, sep="\t", header=T)
seq <- getSeq(hg,nonZerobins$chr,start=nonZerobins$start-flankingLength,end=nonZerobins$end+flankingLength)
seq <- as.character(as.data.frame(seq)[[1]])
nonZerobins.seq <- cbind(nonZerobins,seq)
colnames(nonZerobins.seq) <- c(colnames(nonZerobins),"dna.seq")
fwrite(nonZerobins.seq, file=outputFile, sep="\t", row.names=F, quote=F)
