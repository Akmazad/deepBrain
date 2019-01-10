chrSizeFileName = "hg19.chrom.sizes.txt"
ba9FileName = "normalized_log2_tags_BA9_81_April2015_LR"
ba41FileName = "normalized_log2_tags_BA41_66_Mar2015_LR"
baVermisFileName = "normalized_log2_tags_Vermis_62_Mar2015_LR"
binSize = 200
overlapCutoff = 0.05
bedDir = "/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin"
workingDir = "/Volumes/Data1/PROJECTS/DeepLearning/Test/"
outputFileName = "H3K27ac_binary"

accetylationDat <- function(chrSizeFileName,ba9FileName,ba41FileName,baVermisFileName,binSize,overlapCutoff,bedDir,workingDir,outputFileName){
  setwd(workingDir)
  rm(list=ls())
  
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
  binFile=paste0("hg19_bins_", b,"bp")
  write.table(bins, paste0(binFile,".bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  message("Done",appendLF=T)
  
  ################ Generate bed file of features (H3K27Ac: BA9, BA41, vermis)
  message("Generating bed files for features: ",appendLF=F)
  inDir=workingDir
  outDir=workingDir
  for (feature_file in c(ba9FileName,ba41FileName,baVermisFileName)){
    features=read.csv(paste0(inDir, feature_file, ".csv"))
    #apply(features[,-c(1:3)],2,min)
    features_bed=cbind(features[, c(1:3)], paste(features[,1], features[,2], features[,3], sep="_"), ".")
    write.table(features_bed,paste0(outDir, feature_file, ".bed") , sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  message("Done",appendLF=T)
  
  ############## Overlap Bins with fetures, with a min of 5% overlap ; done in shell using bedTools (can be embeded in R)
  # Step-1: create a shell script namely "AceTylation_Bed_ShellScript.sh" (see attached) within the "workingDir"
  # Step-2: register that script for running with appropriate permission under UNIX using "chmod u+x AceTylation_Bed_ShellScript.sh"
  # Step-3: Put following commands for Bedtools in that shell script which assumes the arguments should be passed from R
  
  ##### "AceTylation_Bed_ShellScript.sh" ########
  ## #!/bin/bash
  ## bedDir=$1	#first argument which is the Bedtools bin directory
  ## cd $2	#second argument which is the working directory
  ## bins=$3	#third argument which is the chromosomal bed file name
  ## overlapCutof=$(echo $4| bc)	#fourth argument which is overlap cut-off
  ## for features in $5 $6 $7	#fifth, sixth, and seventh arguments are b19, b41, and baVermis bedfiles
  ## do
  ## 	  $bedDir/intersectBed -u -f $overlapCutof -a $bins.bed -b $features.bed > $features.overlaps.bed
  ## done
  
  # Step-4: use system2 R function to run this script with arguments passed for the shell script
  message(paste0("Overlapping bins with fetures, with a min of ",overlapCutoff*100, "% overlap: "),appendLF=F)
  system2("./AceTylation_Bed_ShellScript.sh",
            paste(bedDir, 
            workingDir, 
            binFile, 
            overlapCutoff, 
            ba9FileName,
            ba41FileName,
            baVermisFileName,sep=" "))
  # this will create Three overlap bed files
  message("Done",appendLF=T)
    
  ############## Generate the binarised matrix
  message("Generating the binarised matrix: ",appendLF=F)
  setwd(workingDir)
  bins=read.table(paste0(binFile,".bed"), sep="\t", header=FALSE)
  colnames(bins)=c("chr", "start", "end", "id",  "strand")
  feature_files= c(ba9FileName, ba41FileName,baVermisFileName)
  for ( j in c(1:length(feature_files))){
    features=read.csv(paste0(feature_files[j], ".csv"))
    names=colnames(features); rm(features)
    names=names[-c(1:3)]
    overlaps=read.table(paste0(feature_files[j], ".overlaps.bed"))
    colnames(overlaps)=c("chr", "start", "end", "id",  "strand")
    ov=which(bins$id%in%overlaps$id); rm(overlaps)
    binData=matrix(0, nrow=nrow(bins), ncol=length(names))
    colnames(binData)=names
    binData[ov,]=1
    bins=cbind(bins, binData)
    rm(binData)
  }
  write.csv(bins, file=paste0(outputFileName,".csv"), row.names=F)
  message("Done",appendLF=T)
  
}

accetylationDat(chrSizeFileName,ba9FileName,ba41FileName,baVermisFileName,binSize,overlapCutoff,bedDir,workingDir,outputFileName)
