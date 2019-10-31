rm(list = ls(all.names = TRUE))
# Option list
args = commandArgs(trailingOnly=FALSE)
print(args)

dataDir = args[6]
Hu = args[7]
Ep = args[8]
Cg = args[9]
TF = args[10]
seqFile = args[11]
outLabels = args[12]
outCombined = args[13]

setwd(dataDir)
library(data.table)
library(dplyr)

# Example Label data files: 
# 1. HumanFC_nonzero_labels.bed
# 2. EpiMap_nonzero_labels.bed
# 3. CAGE_nonzero_labels.bed
# 4. ENCODE_TFs_nonzero_labels.bed
# 5. HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.bed
# Note, to this point, bins of above files should be the same

dna.dat <- fread(file=seqFile, sep="\t", header=T) # check ids: "grep -o 'e+' HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.dat | wc -l" should return 0
dna.dat$id = paste0(dna.dat$chr, "_", dna.dat$start, "_", dna.dat$end)

human <- fread(Hu, sep="\t", header=T)
human$id <- paste0(human$chr, "_", human$start, "_", human$end) # fix the ids (scientific notation appread!!)
epi <- fread(Ep, sep="\t", header=T)
epi$id <- paste0(epi$chr, "_", epi$start, "_", epi$end)
cg <- fread(Cg, sep="\t", header=T)
cg$id <- paste0(cg$chr, "_", cg$start, "_", cg$end)
tf <- fread(TF, sep="\t", header=T)
tf$id <- paste0(tf$chr, "_", tf$start, "_", tf$end)

output <- cbind(human, epi[which(epi$id %in% human$id),-c(1:4)], cg[which(cg$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
# above line could've just been as follows:
# output <- cbind(human, epi[,-c(1:4)], cg[,-c(1:4)], tf[,-c(1:4)])
output_full <- cbind(dna.dat[which(dna.dat$id %in% human$id), ], human[,-c(1:4)], epi[which(epi$id %in% human$id),-c(1:4)], cg[which(cg$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
colnames(output) <- c(colnames(human), colnames(epi)[-c(1:4)], colnames(cg)[-c(1:4)], colnames(tf)[-c(1:4)])
colnames(output_full) <- c(colnames(dna.dat), colnames(human)[-c(1:4)], colnames(epi)[-c(1:4)], colnames(cg)[-c(1:4)], colnames(tf)[-c(1:4)])
fwrite(output,file=paste0(dataDir,outLabels), sep="\t", row.names=F, quote=F)
fwrite(output_full,file=paste0(dataDir,outCombined), sep="\t", row.names=F, quote=F)
