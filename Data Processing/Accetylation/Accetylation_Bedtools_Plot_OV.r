

accetylationDat <- function(binFile,nBins,chrSizeFileName,ba9FileName,ba41FileName,baVermisFileName,binSize,overlapCutoff,bedDir,workingDir,outputFileName){
  setwd("/Volumes/Data1/PROJECTS/DeepLearning/Test/")
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
  #message("Done",appendLF=T)
  
  #message("Calculating counts: ",appendLF=F)
  setwd(workingDir)
  feature_files= c(ba9FileName, ba41FileName,baVermisFileName)
  max <- -9999
  for (j in c(1:length(feature_files))){
    overlaps=read.table(paste0(feature_files[j], ".overlaps.bed"))
    val <- nrow(overlaps)
    if(val > max) 
      max <- val
  }
  message("Done",appendLF=T)
  retVal <- cbind(max/nBins, 1.0-max/nBins)
  colnames(retVal) <- c("ones","zeros")
  return(retVal)
}  

########################## ONE OFF #################################
chrSizeFileName = "hg19.chrom.sizes.txt"
ba9FileName = "normalized_log2_tags_BA9_81_April2015_LR"
ba41FileName = "normalized_log2_tags_BA41_66_Mar2015_LR"
baVermisFileName = "normalized_log2_tags_Vermis_62_Mar2015_LR"
bedDir = "/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin"
workingDir = "/Volumes/Data1/PROJECTS/DeepLearning/Test/"
outputFileName = "H3K27ac_binary"
binSize = 200
setwd(workingDir)
 
# read in chromosome sizes
chr_size=read.table(chrSizeFileName, sep="\t")
colnames(chr_size)=c("chr", "size")
# remove chromosome patches and sort by chr number
chr_size=chr_size[-grep("_", chr_size$chr, fixed=TRUE),]
chr_size=chr_size[match(paste0("chr", c(c(1:22), "M", "X", "Y")), chr_size$chr), ]

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

################ generate bed file of bins of size b
binSize = 200
workingDir = paste0("/Volumes/Data1/PROJECTS/DeepLearning/Test/",binSize,"bp/")
setwd(workingDir)

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
nBins=nrow(bins)
message("Done",appendLF=T)

  
## iterate over varying overlapCutoff or binSize
dat<-NULL
for(overlapCutoff in seq(0.001,1,0.1)){
  val <- accetylationDat(binFile,nBins,chrSizeFileName,ba9FileName,ba41FileName,baVermisFileName,binSize,overlapCutoff,bedDir,workingDir,outputFileName)
  dat <- rbind(dat,cbind(overlapCutoff,val))
}
write.csv(dat,file=paste0(workingDir,"plot_cutoff_",binSize,"bp.csv"),row.names=F)
