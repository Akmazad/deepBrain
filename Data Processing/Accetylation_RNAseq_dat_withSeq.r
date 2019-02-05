chrSizeFileName = "hg19.chrom.sizes.txt"
ba9FileName = "normalized_log2_tags_BA9_81_April2015_LR"
ba41FileName = "normalized_log2_tags_BA41_66_Mar2015_LR"
baVermisFileName = "normalized_log2_tags_Vermis_62_Mar2015_LR"
rnaSeqFileName = "stringTie.Transcript.SpikeIns_filtered"
binSize = 200
overlapCutoff = 0.05
bedDir = "/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin"
workingDir = "/Volumes/Data1/PROJECTS/DeepLearning/Test/"
outputFileName = "H3K27ac_rnaSeq"
flankingLength=400

Accetylation_RNAseq_dat_withSeq <- function(chrSizeFileName,ba9FileName,ba41FileName,baVermisFileName,rnaSeqFileName,binSize,overlapCutoff,flankingLength,bedDir,workingDir,outputFileName){
  rm(list=ls())
  # get the total human genome
  library("dplyr")
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
  temp_bins=bins[,-4] ## remove the sequence part from Binfile (to save the memory space)
  fwrite(temp_bins, paste0(binFile,".bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  rm(temp_bins); message("Done",appendLF=T)
  
  
  ################ Generate bed file of features (H3K27Ac: BA9, BA41, vermis)
  message("Generating bed files for features: ",appendLF=F)
  inDir=workingDir
  outDir=workingDir
  for (feature_file in c(ba9FileName,ba41FileName,baVermisFileName)){
    features=read.csv(paste0(inDir, feature_file, ".csv"))
    #apply(features[,-c(1:3)],2,min)
    features_bed=cbind(features[, c(1:3)], paste(features[,1], features[,2], features[,3], sep="_"), ".")
   fwrite(features_bed,paste0(outDir, feature_file, ".bed") , sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  message("Done",appendLF=T)
  
  ############## Overlap Bins with fetures, with a min of 5% overlap ; done in shell using bedTools (can be embeded in R)
  # Step-1: create a shell script namely "Ac_RNAseq_Bed_ShellScript.sh" (see attached) within the "workingDir"
  # Step-2: register that script for running with appropriate permission under UNIX using "chmod u+x AceTylation_Bed_ShellScript.sh"
  # Step-3: Put following commands for Bedtools in that shell script which assumes the arguments should be passed from R
  
  ##### "Ac_RNAseq_Bed_ShellScript.sh" ########
  ## #!/bin/bash
  ## bedDir=$1	#first argument which is the Bedtools bin directory
  ## cd $2	#second argument which is the working directory
  ## bins=$3	#third argument which is the chromosomal bed file name
  ## overlapCutof=$(echo $4| bc)	#fourth argument which is overlap cut-off
  ## for features in $5 $6 $7	#fifth, sixth, and seventh arguments are b19, b41, and baVermis bedfiles
  ## do
  ## 	  $bedDir/intersectBed -u -f $overlapCutof -a $bins.bed -b $features.bed > $features.overlaps.bed
  ## done
  ## $bedDir/intersectBed -wa -wb -f $overlapCutof -a $bins.bed -b $8.bed > $8.overlaps.bed

  
  # Step-4: use system2 R function to run this script with arguments passed for the shell script
   message(paste0("Overlapping bins with fetures, with a min of ",overlapCutoff*100, "% overlap: "),appendLF=F)
  system2("./AceTylation_Bed_ShellScript.sh",
            paste(bedDir, 
            workingDir, 
            binFile, 
            overlapCutoff, 
            ba9FileName,
            ba41FileName,
            baVermisFileName,
            sep=" "))
  # this will create Three overlap bed files
  message("Done",appendLF=T)
    
  ############## Generate the binarised matrix
  ##### for Accetylation data
  message("Generating the binarised matrix: ",appendLF=F)
  setwd(workingDir)
  bins=fread(paste0(binFile,".bed"), sep="\t", header=FALSE) ## it is already in the memory
  #colnames(bins)=c("chr", "start", "end", "dna.seq", "id",  "strand")
  feature_files= c(ba9FileName,baVermisFileName,rnaSeqFileName) ## since, ba9 and ba41 are identical
  for ( j in c(1:length(feature_files))){
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
  dir="/Volumes/Data1/PROJECTS/DeepLearning/Test/EncodeDCCExprMatchFiles/"
  feature_files=list.files(dir,pattern="*.narrowPeak$")
  for ( j in c(1:length(feature_files))){
     overlaps=fread(paste0(dir,feature_files[j], ".overlaps.bed"))
     colnames(overlaps)=c("chr", "start", "end", "id",  "strand")
     ov=which(bins$id%in%overlaps$id); rm(overlaps)
     binData=matrix(0, nrow=nrow(bins), ncol=1) ## 1 for each brain region (collapsing all the samples since it's binarized
     colnames(binData)=c(feature_files[j])
     binData[ov,]=1
     bins=cbind(bins, binData)
     rm(binData)
     print(j)
  }

  
  #rnaSeq.features=fread(paste0(rnaSeqFileName, ".csv"),header=T)
  #feature.id=paste(rnaSeq.features$chr,rnaSeq.features$start,rnaSeq.features$end,sep="_")
  #rnaSeq.features=cbind(feature.id,rnaSeq.features[,-c(1:3)])
  #overlaps=fread(paste0(rnaSeqFileName, ".wab.overlaps.bed"))
  #colnames(overlaps)=c("chr", "start", "end", "dna.seq", "bin.id",  "strand", "feature.chr", "feature.start", "feature.end", "feature.id", "feature.strand")
  #overlaps=cbind(overlaps[,overlaps$bin.id],overlaps[,overlaps$feature.id]) ## shortening the overlaps data
  #colnames(overlaps)=c("bin.id","feature.id")
  #overlaps=as.data.frame(overlaps)
  
  ### Do left-outer join teo get 
  #combined <- sort(union(levels(overlaps$feature.id), levels(rnaSeq.features$feature.id)))
  #mx=mutate(overlaps, feature.id=factor(feature.id, levels=combined))
  #my=mutate(rnaSeq.features, feature.id=factor(feature.id, levels=combined))
  #overlaps=left_join(mx,my); rm(combined); rm(mx); rm(my)
  #overlaps=overlaps[,-2] ## get rid of the feature.id, keeping only the bin.id and the sample values
  
  #rnaSeq.features=fread(paste0(rnaSeqFileName, ".overlaps.bed"),header=F)
  #rnaSeq.features=rnaSeq.features[,-c(4:11)]
  ## union of similar row
  
  ## find overlap bins
  #ov=which(bins$id%in%overlaps$bin.id); rm(overlaps)
  
  ## binarize
  
  ## rbind with the main output 
  
  #ov=which(bins$id%in%overlaps$bin.id); rm(overlaps)
  #binData=matrix(0, nrow=nrow(bins), ncol=1)
  
  #write.csv(bins, file=paste0(outputFileName,".csv"), row.names=F)
  fwrite(bins, file=paste0(outputFileName,".csv"), row.names=F, quote=F)
  message("Done",appendLF=T)

}

Accetylation_RNAseq_dat_withSeq(chrSizeFileName,ba9FileName,ba41FileName,baVermisFileName,rnaSeqFileName,binSize,overlapCutoff,flankingLength,bedDir,workingDir,outputFileName)
