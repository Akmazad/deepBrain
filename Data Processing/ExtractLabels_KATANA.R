rm(list = ls(all.names = TRUE))
# Option list
args = commandArgs(trailingOnly=FALSE)
print(args)

dataDir = args[6]
H = args[7]
E = args[8]
T = args[9]
seqFile = args[10]
outLabels = args[11]
outCombined = args[12]

setwd(dataDir)
library(data.table)
library(dplyr)

# Example Label data files: 
# 1. HumanFC_nonzero_labels.bed
# 2. EpiMap_nonzero_labels.bed
# 3. ENCODE_TFs_nonzero_labels.bed
# 4. HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.bed

dna.dat <- fread(file=seqFile, sep="\t", header=T) # ids seemed fine: "grep -o 'e+' HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.dat | wc -l" return 0

human <- fread(H, sep="\t", header=T)
human$id <- paste0(human$chr, "_", human$start, "_", human$end) # fix the ids (scientific notation appread!!)
epi <- fread(E, sep="\t", header=T)
epi$id <- paste0(epi$chr, "_", epi$start, "_", epi$end) # fix the ids (scientific notation appread!!)
tf <- fread(T, sep="\t", header=T)
tf$id <- paste0(tf$chr, "_", tf$start, "_", tf$end) # fix the ids (scientific notation appread!!)

output <- cbind(human, epi[which(epi$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
output_full <- cbind(dna.dat[which(dna.dat$id %in% human$id), ], human[,-c(1:4)], epi[which(epi$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
colnames(output) <- c(colnames(human), colnames(epi)[-c(1:4)], colnames(tf)[-c(1:4)])
colnames(output_full) <- c(colnames(dna.dat), colnames(human)[-c(1:4)], colnames(epi)[-c(1:4)], colnames(tf)[-c(1:4)])
fwrite(output,file=outLabels, sep="\t", row.names=F, quote=F)
fwrite(output_full,file=outCombined, sep="\t", row.names=F, quote=F)