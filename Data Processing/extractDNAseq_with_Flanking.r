chrSizeFileName = "hg19.chrom.sizes.txt"
binSize = 200
flankingLength=400
workingDir = "/Volumes/Data1/PROJECTS/DeepLearning/Test/"

Prep_Acc_RNAseq_dat_withSeq <- function(chrSizeFileName,binSize,flankingLength,workingDir){
  # get the total human genome
  library(dplyr)
  library("data.table")
  library("BSgenome.Hsapiens.UCSC.hg19")
  hg <- BSgenome.Hsapiens.UCSC.hg19
  
  setwd(workingDir)
  b=binSize
  binFile=paste0("hg19_bins_", b,"bp")
  
  # read in chromosome sizes
  chr_size=read.table(chrSizeFileName, sep="\t")
  colnames(chr_size)=c("chr", "size")
  # remove chromosome patches and sort by chr number
  chr_size=chr_size[-grep("_", chr_size$chr, fixed=TRUE),]
  chr_size=chr_size[match(paste0("chr", c(c(1:22), "M", "X", "Y")), chr_size$chr), ]
  
  # generate bed file of bins of size b
  message("Generating bed files for each bins of size b: ",appendLF=F)
  for (j in c(1:nrow(chr_size))){
    # these start and end positions doesn't consider initial few bins that don't flanking sequence
    start=seq(from=flankingLength, to=chr_size$size[j]-flankingLength-b, by=b)+1 
    end=seq(from=flankingLength+b, to=chr_size$size[j]-flankingLength, by=b) 
    
    # retrieve dna subsequence (with flanking) from hg
    chrName=as.character(chr_size$chr[j])
    dnaSeq.start=start-flankingLength
    dnaSeq.end=end+flankingLength
    fasta.seq=getSeq(hg,chrName,start=dnaSeq.start,end=dnaSeq.end)
    tempFasta = as.character(as.data.frame(fasta.seq)[[1]])
    
    chr_bins=cbind(chrName,start[1:length(end)],end)
    chr_bins=cbind(chr_bins,tempFasta)
    if (j==1) bins=chr_bins else bins=rbind(bins, chr_bins) 
    print(j)
  }
  bins=as.data.frame(bins)
  colnames(bins)=c("chr", "start", "end", "dna.seq")
  bins$id=paste(bins$chr, bins$start, bins$end, sep="_")
  
  fwrite(bins,file=paste0(binFile,".fasta.bed"), sep="\t", row.names=F, quote=F)

}
