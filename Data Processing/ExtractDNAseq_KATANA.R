rm(list = ls(all.names = TRUE))
# Option list
args = commandArgs(trailingOnly=FALSE)
print(args)

dataDir = args[6]
flankingLength = strtoi(args[7])
inputFile = args[8]
outputFile = args[9]

setwd(dataDir)
library(data.table)
library(dplyr)

library("BSgenome.Hsapiens.UCSC.hg19")
hg <- BSgenome.Hsapiens.UCSC.hg19

# read non-zero bins
nonZerobins <- fread(inputFile, sep="\t", header=T)
if(prod(colnames(nonZerobins) == c("chr", "start", "end", "binID")) == 0){
	colnames(nonZerobins) <- c("chr", "start", "end", "binID")
}

seq <- getSeq(hg,nonZerobins$chr,start=nonZerobins$start-flankingLength,end=nonZerobins$end+flankingLength)
seq <- as.character(as.data.frame(seq)[[1]])
nonZerobins.seq <- cbind(nonZerobins,seq)
colnames(nonZerobins.seq) <- c(colnames(nonZerobins),"dna.seq")
fwrite(nonZerobins.seq, file=outputFile, sep="\t", row.names=F, quote=F)
