# Downloaded hg19 chromosome sizes from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
setwd("/Volumes/Data1/PROJECTS/DeepLearning/Test/")
rm(list=ls())
# read in chromosome sizes
chr_size=read.table("hg19.chrom.sizes.txt", sep="\t")
colnames(chr_size)=c("chr", "size")
# remove chromosome patches and sort by chr number
chr_size=chr_size[-grep("_", chr_size$chr, fixed=TRUE),]
chr_size=chr_size[match(paste0("chr", c(c(1:22), "M", "X", "Y")), chr_size$chr), ]
################ generate bed file of bins of size b
b=200
for (j in c(1:nrow(chr_size)))
{
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
write.table(bins, paste0("hg19_bins_", b,"bp.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

################ Generate bed file of features (H3K27Ac: BA9, BA41, vermis)
inDir="/Volumes/Data1/PROJECTS/DeepLearning/Test/"
outDir="/Volumes/Data1/PROJECTS/DeepLearning/Test/"
for (feature_file in c("normalized_log2_tags_BA9_81_April2015_LR", "normalized_log2_tags_BA41_66_Mar2015_LR","normalized_log2_tags_Vermis_62_Mar2015_LR"))
{
features=read.csv(paste0(inDir, feature_file, ".csv"))
#apply(features[,-c(1:3)],2,min)
features_bed=cbind(features[, c(1:3)], paste(features[,1], features[,2], features[,3], sep="_"), ".")
write.table(features_bed,paste0(outDir, feature_file, ".bed") , sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

############## Overlap Bins with fetures, with a min of 5% overlap ; done in shell using bedTools (can be embeded in R)
runDir=/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin
cd /Volumes/Data1/PROJECTS/DeepLearning/Test
bins=hg19_bins_200bp
for features in "normalized_log2_tags_BA9_81_April2015_LR" "normalized_log2_tags_BA41_66_Mar2015_LR" "normalized_log2_tags_Vermis_62_Mar2015_LR"
do
$runDir/intersectBed -u -f 0.05 -a $bins.bed  -b $features.bed > $features.overlaps.bed
done

############## Generate the binarised matrix
rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/DeepLearning/Test")

bins=read.table("hg19_bins_200bp.bed", sep="\t", header=FALSE)
colnames(bins)=c("chr", "start", "end", "id",  "strand")
feature_files= c("normalized_log2_tags_BA9_81_April2015_LR", "normalized_log2_tags_BA41_66_Mar2015_LR","normalized_log2_tags_Vermis_62_Mar2015_LR")
for ( j in c(1:length(feature_files)))
{
  features=read.csv(paste0(feature_files[j], ".csv"))
  names=colnames(features); rm(features)
  names=names[-c(1:3)]
  overlaps=read.table(paste0(feature_files[j], ".bed"))
  colnames(overlaps)=c("chr", "start", "end", "id",  "strand")
  ov=which(bins$id%in%overlaps$id); rm(overlaps)
  binData=matrix(0, nrow=nrow(bins), ncol=length(names))
  colnames(binData)=names
  binData[ov,]=1
  bins=cbind(bins, binData)
  rm(binData)
}
save(bins, file="H3K27ac_binary.rda")