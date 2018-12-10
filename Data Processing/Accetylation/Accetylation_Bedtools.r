# chrSizeFileName = 
# AcetylationFileName = 
# binSize = 
# overlapCutoff = 
# bedDir = 
# workingDir = 

accetylationDat <- function(chrSizeFileName,AcetylationFileName,binSize,overlapCutoff,bedDir,workingDir){
  setwd(workingDir)
  rm(list=ls())
  
  # read in chromosome sizes
  chr_size=read.table(chrSizeFileName, sep="\t")
  colnames(chr_size)=c("chr", "size")
  # remove chromosome patches and sort by chr number
  chr_size=chr_size[-grep("_", chr_size$chr, fixed=TRUE),]
  chr_size=chr_size[match(paste0("chr", c(c(1:22), "M", "X", "Y")), chr_size$chr), ]

  ################ generate bed file of bins of size b
  b=binSize
  for (j in c(1:nrow(chr_size))){
    start=seq(from=0, to=chr_size$size[j], by=b)+1
    end=seq(from=b, to=chr_size$size[j], by=b)
    chr_bins=cbind(as.character(chr_size$chr[j]),start[1:length(end)],end)
    if (j==1) bins=chr_bins else bins=rbind(bins, chr_bins) 
    print(j)
  }
  bins=as.data.frame(bins)
  colnames(bins)=c("chr", "start", "end")
  bins$id=paste(bins$chr, bins$start, bins$end, sep="_")
  bins$strand="."
  binFile=paste0("hg19_bins_", b,"bp.bed")
  write.table(bins, binFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

}
