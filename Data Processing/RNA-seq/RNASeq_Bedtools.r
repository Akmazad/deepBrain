chrSizeFileName = "hg19.chrom.sizes.txt"
rnaSeqFileName = "StringTie_humanTFsOnly"
binSize = 200
overlapCutoff = 0.05
bedDir = "/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin"
workingDir = "/Volumes/Data1/PROJECTS/DeepLearning/StringTie/Files/"
outputFileName = "rnaSeq_binary"

accetylationDat <- function(chrSizeFileName,rnaSeqFileName,binSize,overlapCutoff,bedDir,workingDir,outputFileName){
  setwd(workingDir)
  #rm(list=ls())
  
  # read in chromosome sizes
  chr_size=read.table(chrSizeFileName, sep="\t")
  colnames(chr_size)=c("chr", "size")
  # remove chromosome patches and sort by chr number
  chr_size=chr_size[-grep("_", chr_size$chr, fixed=TRUE),]
  chr_size=chr_size[match(paste0("chr", c(c(1:22), "M", "X", "Y")), chr_size$chr), ]

  ################ generate bed file of bins of size b
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
  bins$strand="."
  binFile=paste0("hg19_bins_", b,"bp.bed")
  write.table(bins, binFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  message("Done",appendLF=T)
  
  ################ Generate bed file of features from RNAseq
  message("Generating bed files for features: ",appendLF=F)
  inDir=workingDir
  outDir=workingDir
  feature_file=rnaSeqFileName
  features=read.csv(paste0(inDir, rnaSeqFileName, ".csv"))
  features[,3] = paste0("chr",features[,3])
  features_bed=cbind(features[,3],features[,5],features[,6], paste(features[,3], features[,5], features[,6], sep="_"), ".")
  write.table(features_bed,paste0(outDir, feature_file, ".bed") , sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  message("Done",appendLF=T)
  
  ############## Overlap Bins with fetures, with a min of 5% overlap ; done in shell using bedTools (can be embeded in R)
  # Use 'system2' R function to run this bedtool command with arguments passed for the shell script
  message(paste0("Overlapping bins with fetures, with a min of ",overlapCutoff*100, "% overlap: "),appendLF=F)
  system2("/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin/intersectBed",
          paste("-u -f", 
            overlapCutoff, 
            "-a",
            paste0(outDir, binFile),
            "-b",
            paste0(outDir, feature_file, ".bed"),
            ">",
            paste0(outDir, feature_file, ".overlaps.bed"),
            sep=" "))
  # this will create Three overlap bed files
  message("Done",appendLF=T)
  
  ############## Generate the binarised matrix
  message("Generating the binarised matrix: ",appendLF=F)
  setwd(workingDir)
  bins=read.table(binFile, sep="\t", header=FALSE)
  colnames(bins)=c("chr", "start", "end", "id",  "strand")
  features=read.csv(paste0(feature_file, ".csv"))
  names=colnames(features); rm(features)
  names=names[-c(1:11)]
  overlaps=read.table(paste0(feature_file, ".overlaps.bed"))
  colnames(overlaps)=c("chr", "start", "end", "id",  "strand")
  ov=which(bins$id%in%overlaps$id); rm(overlaps)
  binData=matrix(0, nrow=nrow(bins), ncol=length(names))
  colnames(binData)=names
  binData[ov,]=1
  bins=cbind(bins, binData)
  rm(binData)
}
