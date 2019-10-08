rm(list = ls(all.names = TRUE))
# Option list
args = commandArgs(trailingOnly=FALSE)
print(args)

dataDir = args[6]
Hu = args[7]
Ep = args[8]
TF = args[9]
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

dna.dat <- fread(file=seqFile, sep="\t", header=T) # check ids: "grep -o 'e+' HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.dat | wc -l" should return 0
dna.dat$id = paste0(dna.dat$chr, "_", dna.dat$start, "_", dna.dat$end)

human <- fread(Hu, sep="\t", header=T)
human$id <- paste0(human$chr, "_", human$start, "_", human$end) # fix the ids (scientific notation appread!!)
epi <- fread(Ep, sep="\t", header=T)
epi$id <- paste0(epi$chr, "_", epi$start, "_", epi$end)
tf <- fread(TF, sep="\t", header=T)
tf$id <- paste0(tf$chr, "_", tf$start, "_", tf$end)

output <- cbind(human, epi[which(epi$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
output_full <- cbind(dna.dat[which(dna.dat$id %in% human$id), ], human[,-c(1:4)], epi[which(epi$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
colnames(output) <- c(colnames(human), colnames(epi)[-c(1:4)], colnames(tf)[-c(1:4)])
colnames(output_full) <- c(colnames(dna.dat), colnames(human)[-c(1:4)], colnames(epi)[-c(1:4)], colnames(tf)[-c(1:4)])
fwrite(output,file=paste0(dataDir,outLabels), sep="\t", row.names=F, quote=F)
fwrite(output_full,file=paste0(dataDir,outCombined), sep="\t", row.names=F, quote=F)
