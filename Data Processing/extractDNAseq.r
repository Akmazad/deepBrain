chrSizeFileName = "hg19.chrom.sizes.txt"
workingDir = "/Volumes/Data1/PROJECTS/DeepLearning/Test/"
library("BSgenome.Hsapiens.UCSC.hg19")
hg <- BSgenome.Hsapiens.UCSC.hg19
setwd(workingDir)

extractDNAseq <- function(hg,chrSizeFileName,binSize,flankingLength,workingDir){
  # read in chromosome sizes
  chr_size=read.table(chrSizeFileName, sep="\t")
  colnames(chr_size)=c("chr", "size")
  # remove chromosome patches and sort by chr number
  chr_size=chr_size[-grep("_", chr_size$chr, fixed=TRUE),]
  chr_size=chr_size[match(paste0("chr", c(c(1:22), "M", "X", "Y")), chr_size$chr), ]

  ################ generate bed file of bins of size b
  message(paste0("Generating bed files for each bins of size b: ",binSize),appendLF=F)
  b=binSize
  for (j in c(1:nrow(chr_size))){
    start=seq(from=0, to=chr_size$size[j]-b, by=b)+1 - flankingLength
    end=seq(from=b, to=chr_size$size[j], by=b) + flankingLength
    chrName=as.character(chr_size$chr[j])
    fasta.seq=getSeq(hg,chrName,start=start,end=end)
    tempFasta = as.character(as.data.frame(fasta.seq)[[1]])
    chr_bins=cbind(chrName,start[1:length(end)],end)
    chr_bins=cbind(chr_bins,tempFasta)
    if (j==1) bins=chr_bins else bins=rbind(bins, chr_bins) 
  }
  bins=as.data.frame(bins)
  colnames(bins)=c("chr", "start", "end","FastaSeq")
  bins$id=paste(bins$chr, bins$start, bins$end, sep="_")
  bins$strand="."
  binFile=paste0("hg19.binwise.fasta.", b,"bp")
  write.table(bins, paste0(binFile,".bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  message("Done",appendLF=T)
}

## run this code for a set of bins
flakingLength=1000
for(binSize in seq(200,6400,by=200)){
  extractDNAseq(hg,chrSizeFileName,binSize,flakingLength,workingDir)
}
## for try
#flakingLength=1000
#binSize = 200
#j=1
#b=binSize
#start=seq(from=0, to=chr_size$size[j]-b, by=b)+1
#end=seq(from=b, to=chr_size$size[j], by=b)
#chrName=as.character(chr_size$chr[j])
#fasta.seq=getSeq(hg,chrName,start=start,end=end)
#tempFasta = as.character(as.data.frame(fasta.seq)[[1]])
#chr_bins=cbind(chrName,start[1:length(end)],end)
#chr_bins=cbind(chr_bins,tempFasta)
#colnames(chr_bins)=c("chr", "start", "end","FastaSeq")
#bins=as.data.frame(bins)
#colnames(bins)=c("chr", "start", "end","FastaSeq")
#binFile=paste0("hg19.binwise.fasta.", b,"bp")
#write.table(bins, paste0(binFile,".bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
