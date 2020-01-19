# Input:  binsize, flankingLength, howManyBins (default: 1M)
#         mergedPeaks (mergedPeakHeightMatrix_HumanFC_filtered.bed, 
#         mergedPeakHeightMatrix_EpiMap_filtered.bed, 
#         final.dat.tf.overlaps.bed)
# Host: RNA machine/KATANA 
#       RNA:    Data dir: /Volumes/Data1/PROJECTS/DeepLearning/Test
#               Script dir: /Volumes/Data1/PROJECTS/DeepLearning/Test
#       KATANA: Data dir: /srv/scratch/z3526914/DeepBrain/Data
#               Script dir: /srv/scratch/z3526914/DeepBrain/Scripts
options(scipen=999) # prevent auto-scientific notations of numeric values

library(dplyr)
library(data.table)

# ------------- Set Bin-flanking configuration
binSize = 200
flanking = 4*binSize     # can be arbitrarily given
dataDir = paste0("/Volumes/Data1/PROJECTS/DeepLearning/Test/",binSize,"_",flanking,"/")
# dataDir = paste0("/srv/scratch/z3526914/DeepBrain/Data/",binSize,"_",flanking,"/")

dir.create(dataDir)
# set a fixed number of bins for this bin_flanking test:
# Rationale: for smaller bin/flanking size, the number of bins will be huge
#            for which downstream data processing may suffer resource issue
howManyBins = 1000000
# chrSizeFileName = "/srv/scratch/z3526914/DeepBrain/Data/hg19.chrom.sizes.txt"
chrSizeFileName = "/Volumes/Data1/PROJECTS/DeepLearning/Test/hg19.chrom.sizes.txt"
chr_size=read.table(chrSizeFileName, sep="\t")
colnames(chr_size)=c("chr", "size")
# remove chromosome patches and sort by chr number
chr_size=chr_size[-grep("_", chr_size$chr, fixed=TRUE),]
chr_size=chr_size[match(paste0("chr", c(c(1:22), "M", "X", "Y")), chr_size$chr), ]

# 1. generate binIDs with size given, and chose randomly 1M of them: one File output
# generate bed file of bins of size b
message("Generating bed files for each bins of size b: ",appendLF=F)
b=binSize
for (j in c(1:nrow(chr_size))){
  start=seq(from=0, to=chr_size$size[j], by=b)+1
  end=seq(from=b, to=chr_size$size[j], by=b)
  chr_bins=cbind(as.character(chr_size$chr[j]),start[1:length(end)],end)
  if (j==1) bins=chr_bins else bins=rbind(bins, chr_bins) 
}
bins=as.data.frame(bins)
colnames(bins)=c("chr", "start", "end")
bins$id=paste(bins$chr, bins$start, bins$end, sep="_")
binFile=paste0("hg19_bins_", b,"bp")
fwrite(bins, file=paste0(dataDir,binFile,".bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
message("Done",appendLF=T)
# select randomly a fixed number of bins
bin_inputFile = binFile
bin_outputFile = paste0(bin_inputFile,"_rand")

message(paste("Select", howManyBins, "bins randomly "),appendLF=F)
system2('shuf', 
        paste('-n', howManyBins, paste0(dataDir,bin_inputFile,".bed"), '>',paste0(dataDir,bin_outputFile,".bed"), sep=' '), 
        wait=T)
message("Done",appendLF=T)

