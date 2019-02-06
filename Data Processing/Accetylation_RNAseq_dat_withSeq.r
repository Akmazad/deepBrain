chrSizeFileName = "hg19.chrom.sizes.txt"
ba9FileName = "normalized_log2_tags_BA9_81_April2015_LR"
ba41FileName = "normalized_log2_tags_BA41_66_Mar2015_LR"
baVermisFileName = "normalized_log2_tags_Vermis_62_Mar2015_LR"
rnaSeqFileName = "stringTie.Transcript.SpikeIns_filtered"
binSize = 200
overlapCutoff = 0.05
flankingLength=400

# for "rna" machine
#bedDir = "/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin"
#workingDir = "/Volumes/Data1/PROJECTS/DeepLearning/Test/"
#outputFileName = "H3K27ac_rnaSeq"

# for Raijin
workingDir = "/short/yr31/aa7970/azData/DeepBrain/Data/"
Pos.OutputFileName = "H3K27ac_rnaSeq.Pos.dat"
Neg.OutputFileName = "H3K27ac_rnaSeq.Neg.dat"
Comb.OutputFileName= "H3K27ac_rnaSeq.Pos.Neg.dat"


Prep_Acc_RNAseq_dat_withSeq <- function(chrSizeFileName,ba9FileName,ba41FileName,baVermisFileName,rnaSeqFileName,binSize,overlapCutoff,flankingLength,workingDir,outputFileName){
  rm(list=ls())
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
  
  ################ generate bed file of bins of size b
  message("Generating bed files for each bins of size b: ",appendLF=F)
  for (j in c(1:nrow(chr_size))){
    ### these start and end positions doesn't consider initial few bins that don't flanking sequence
    start=seq(from=flankingLength, to=chr_size$size[j]-flankingLength-b, by=b)+1 
    end=seq(from=flankingLength+b, to=chr_size$size[j]-flankingLength, by=b) 
    
    ## retrieve dna subsequence (with flanking) from hg
    chrName=as.character(chr_size$chr[j])
    dnaSeq.start=start-flankingLength
    dnaSeq.end=end+flankingLength
    fasta.seq=getSeq(hg,chrName,start=dnaSeq.start,end=dnaSeq.end)
    tempFasta = as.character(as.data.frame(fasta.seq)[[1]])
    
    chr_bins=cbind(chrName,start[1:length(end)],end)
    chr_bins=cbind(chr_bins,tempFasta)
    if (j==1) bins=chr_bins else bins=rbind(bins, chr_bins) 
  }
  bins=as.data.frame(bins)
  colnames(bins)=c("chr", "start", "end", "dna.seq")
  bins$id=paste(bins$chr, bins$start, bins$end, sep="_")
  bins$strand="+" ## by default the strand is "+"
  #temp_bins=bins[,-4] ## remove the sequence part from Binfile (to save the memory space)
  #fwrite(temp_bins, paste0(binFile,".bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  #rm(temp_bins); 
  message("Done",appendLF=T)
  
  
  ############## Generate the binarised matrix
  ##### for Accetylation data
  message("Generating the binarised matrix: ",appendLF=F)
  setwd(workingDir)
  #bins=fread(paste0(binFile,".bed"), sep="\t", header=FALSE) ## it is already in the memory
  #colnames(bins)=c("chr", "start", "end", "dna.seq", "id",  "strand")
  feature_files= c(ba9FileName,baVermisFileName,rnaSeqFileName) ## since, ba9 and ba41 are identical
  for (j in c(1:length(feature_files))){
    overlaps=fread(paste0(feature_files[j], ".overlaps.bed"))
    colnames(overlaps)=c("chr", "start", "end", "id",  "strand")
    ov=which(bins$id%in%overlaps$id); rm(overlaps)
    binData=matrix(0, nrow=nrow(bins), ncol=1) ## 1 for each brain region (collapsing all the samples since it's binarized
    colnames(binData)=c(feature_files[j])
    binData[ov,]=1
    bins=cbind(bins, binData)
    rm(binData)
  }
  
  ##### for RNA-seq data
  dir=paste0(workingDir, "tfbs_intersectBED_overlap_files/")
  feature_files=list.files(dir,pattern="*.narrowPeak.overlaps.bed$")
  for (j in c(1:length(feature_files))){
    overlaps=fread(paste0(dir,feature_files[j]))
    colnames(overlaps)=c("chr", "start", "end", "id",  "strand")
    ov=which(bins$id%in%overlaps$id); rm(overlaps)
    binData=matrix(0, nrow=nrow(bins), ncol=1) ## 1 for each brain region (collapsing all the samples since it's binarized
    colnames(binData)=c(feature_files[j])
    binData[ov,]=1
    bins=cbind(bins, binData)
    rm(binData)
    print(j)
  }
  
  ind=7
  pos.ind=which(apply(bins, 1, FUN = function(x) any(x[c(ind:ncol(bins))]==1)))
  pos.bins = bins[pos.ind,]
  neg.bins = bins[-pos.ind,]
  pos.bins = if(nrow(pos.bins) > nrow(neg.bins)) pos.bins[sample(nrow(pos.bins),replace = F, nrow(neg.bins)),] else pos.bins
  neg.bins = if(nrow(neg.bins) > nrow(pos.bins)) neg.bins[sample(nrow(neg.bins),replace = F, nrow(pos.bins)),] else neg.bins
  
  fwrite(pos.bins, file=Pos.OutputFileName, row.names=F, quote=F, sep="\t")
  fwrite(neg.bins, file=Neg.OutputFileName, row.names=F, quote=F, sep="\t")
  fwrite(rbind(pos.bins,neg.bins), file=Comb.OutputFileName, row.names=F, quote=F, sep="\t")
                      
  message("Done",appendLF=T)
  
}

Prep_Acc_RNAseq_dat_withSeq(chrSizeFileName,ba9FileName,ba41FileName,baVermisFileName,rnaSeqFileName,binSize,overlapCutoff,flankingLength,workingDir,outputFileName)
