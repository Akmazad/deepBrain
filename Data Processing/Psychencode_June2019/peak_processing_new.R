### Author: Irina Voineagu, Assoc Professor, School of BABS, UNSW Sydney
### Modified: AKM Azad, Postdoc, School of BABS, UNSW Sydeny


#######################STEP1: For each sample, add a column with the sample name and one with peakCoordinate_sample
rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/Psychencode_June2019/BrainGVEX/ATACseq_Peaks/")
# read in list of peak fiels and metadata
files=read.table("fileList.txt", sep="\t", header=FALSE)
meta=read.table("../Metadata/SYNAPSE_METADATA_MANIFEST.tsv",sep="\t", header=TRUE)
# runtime for 292 samples: 15min
for (j in c(1:nrow(files)))
{
  data=read.table(as.character(files[j,1]), sep="\t", header=FALSE)
  sample=meta$specimenID[grep (files[j,1], meta$path, fixed = T)]
  data=cbind(data, sample, paste(data[,1],data[,2],data[,3], sample, sep="_"))
  # data=cbind(data, sample, paste(data[,1],data[,2],data[,3], sample, sep="&"))    # For EpiMap
  write.table(data, paste0(files[j,1],".sampleID.bed"), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
  print(j)
}

#######################STEP2: Merge peaks using BEDTOOLS
cd /Volumes/Data1/PROJECTS/Psychencode_June2019/BrainGVEX/ATACseq_Peaks/
  outdir=/Volumes/Data1/PROJECTS/Psychencode_June2019/BrainGVEX/ProcessedData/
  # concatenate all peak files in one; runtime 3min
  cat *.sampleID.bed >> $outdir/allPeaks.sampleID.bed
#wc -l $outdir/allPeaks.sampleID.bed
# 11687419 
# sort and merge peaks (mergebed requires sorted input)
# sort: runtime:4 min
sortbed -i $outdir/allPeaks.sampleID.bed > $outdir/allPeaks.sampleID.sorted.bed
# merge: runtime: 3min
mergebed -c 12 -o collapse -i $outdir/allPeaks.sampleID.sorted.bed > $outdir/allPeaks.sampleID.merged.bed
# mergebed -c 11 -o collapse -i $outdir/allPeaks.sampleID.sorted.bed > $outdir/allPeaks.sampleID.merged.bed     # for Yale_ASD
# wc -l $outdir/allPeaks.sampleID.merged.bed
#1049296
#######################STEP3: Generate a data matrix with peak height info for merged peaks; Runtime ~ 30min
library(tidyr)
library(dplyr)
library(data.table)
rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/Psychencode_June2019/BrainGVEX/ATACseq_Peaks/")
#read in merged peak coordinates
merged=read.table("/Volumes/Data1/PROJECTS/Psychencode_June2019/BrainGVEX/ProcessedData/allPeaks.sampleID.merged.bed", sep="\t")
# read in list of peak files and metadata
files=read.table("/Volumes/Data1/PROJECTS/Psychencode_June2019/BrainGVEX/ATACseq_Peaks/fileList.txt", sep="\t", header=FALSE)
meta=read.table("/Volumes/Data1/PROJECTS/Psychencode_June2019/BrainGVEX/Metadata/SYNAPSE_METADATA_MANIFEST.tsv",sep="\t", header=TRUE)

# uncollapse merged peak table to one entry per row
merged$V4=as.character(merged$V4)
merged_uncollapsed= merged %>% unnest(V4 = strsplit(V4, ","))
dim(merged_uncollapsed)
#[1] 11687419        4

# remove duplicate rows. not sure why they exist
merged_uncollapsed=distinct(merged_uncollapsed)
dim(merged_uncollapsed)
#[1] 9709119 

# split column 4 into initialPeakID and sample name and add a mergedPeakID
initialPeakInfo=transpose(strsplit(merged_uncollapsed$V4, "_"))
# initialPeakInfo=transpose(strsplit(merged_uncollapsed$V4, "&"))    # For EpiMap
merged_uncollapsed=cbind(merged_uncollapsed, 
                         paste(merged_uncollapsed$V1, merged_uncollapsed$V2, merged_uncollapsed$V3, sep="_"),
                         paste(initialPeakInfo[[1]], initialPeakInfo[[2]], initialPeakInfo[[3]], sep="_"),
                         initialPeakInfo[[4]])

colnames(merged_uncollapsed)=c("chr", "start", "end", "initialPeakID_Sample", "mergedPeakID", "initialPeakID", "Sample")
# make a placeholder matrix
peaks=matrix(NA, nrow=length(unique(merged_uncollapsed$mergedPeakID)), ncol=nrow(meta))
rownames(peaks)=unique(merged_uncollapsed$mergedPeakID)
colnames(peaks)=as.character(meta$specimenID)
# specify peakHeight column from the initial peak files (can differ between datasets), as well as initialPeak_SampleCol and sampleCol  (which were created in step 1).
peakHeightCol=7
sampleCol=11  # 10 for Yale_ASD
initialPeak_SampleCol=12  # 11 for Yale_ASD
# for each sample
for (j in c(1:nrow(files)))
{ # read in peak height data
  data=read.table(paste0(as.character(files[j,1]), ".sampleID.bed"), sep="\t", header=FALSE)
  # subset merged peaks by those with contributing peaks from sample j
  use_peaks=merged_uncollapsed[which(merged_uncollapsed$Sample%in%data[,sampleCol]),]
  # add peak height info
  use_peaks$peakHeight=data[match(use_peaks$initialPeakID_Sample, data[,initialPeak_SampleCol]), peakHeightCol]
  # sum all peaks from sampleJ by merged peak coordinates
  peak_sum=aggregate(use_peaks$peakHeigh, by=as.data.frame(use_peaks$mergedPeakID), FUN="sum")
  # add data to the matrix            
  s=which(colnames(peaks)%in%data[, sampleCol])
  p=match(rownames(peaks), peak_sum[,1])
  peaks[,s]=peak_sum[p,2]
  print(j)
}
save(peaks,file= "/Volumes/Data1/PROJECTS/Psychencode_June2019/BrainGVEX/ProcessedData/mergedPeakHeightMatrix.rda")
