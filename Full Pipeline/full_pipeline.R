# Input:  binsize, flankingLength, howManyBins (default: 1M)
#         mergedPeaks (mergedPeakHeightMatrix_HumanFC_filtered.bed, 
#         mergedPeakHeightMatrix_EpiMap_filtered.bed, 
#         Brain_CagePeaks_filtered.bed,
#         final.tf.overlaps.bed)
# Host: RNA machine/KATANA 
#       RNA:    Data dir: /Volumes/Data1/PROJECTS/DeepLearning/Test
#               Script dir: /Volumes/Data1/PROJECTS/DeepLearning/Test
#       KATANA: Data dir: /srv/scratch/z3526914/DeepBrain/Data
#               Script dir: /srv/scratch/z3526914/DeepBrain/Scripts
options(scipen=999) # prevent auto-scientific notations of numeric values

binSize = 200
flanking = 4*binSize     # can be arbitrarily given
dataDir = paste0("/Volumes/Data1/PROJECTS/DeepLearning/Test/",binSize,"_",flanking,"/")
# dataDir = paste0("/srv/scratch/z3526914/DeepBrain/Data/",binSize,"_",flanking,"/")

dir.create(dataDir, recursive=T)
# set a fixed number of bins for this bin_flanking test:
# Rationale: for smaller bin/flanking size, the number of bins will be huge
#            for which downstream data processing may suffer resource issue
howManyBins = 1000000
# chrSizeFileName = "/srv/scratch/z3526914/DeepBrain/Data/hg19.chrom.sizes.txt"
chrSizeFileName = "/Volumes/Data1/PROJECTS/DeepLearning/Test/hg19.chrom.sizes.txt"
chr_size=read.table(chrSizeFileName, sep="\t")
colnames(chr_size)=c("chr", "size")
# remove chromosome patches and sort by chr number
chr_size=chr_size[-grep("_", chr_size$chr, fixed=TRUE),]
chr_size=chr_size[match(paste0("chr", c(c(1:22), "M", "X", "Y")), chr_size$chr), ]

library(dplyr)
library(data.table)

# 1. generate binIDs with size given, and chose randomly 1M of them: one File output
# generate bed file of bins of size b
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
# Note: no strand is mentioned, hence 1 less column in all the subsequent files

binFile=paste0("hg19_bins_", b,"bp")
fwrite(bins, file=paste0(dataDir,binFile,".bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
message("Done",appendLF=T)
# select randomly a fixed number of bins
bin_inputFile = binFile
bin_outputFile = paste0(bin_inputFile,"_rand")

message(paste("Select", howManyBins, "bins randomly "),appendLF=F)
system2('shuf', 
        paste('-n', howManyBins, paste0(dataDir,bin_inputFile,".bed"), '>',paste0(dataDir,bin_outputFile,".bed"), sep=' '), 
        wait=T)
message("Done",appendLF=T)

# 2. IntersectBED bins (-a) with mergedPeaks (-b) in shell and post-processing: three file outputs

# ---------------------------------------------- HumanFC
humanFC_inputFile = "mergedPeakHeightMatrix_HumanFC_filtered"
# Note: "humanFC_inputFile" should exist within the "dataDir", and others (Epimap, Cage, and TFs) are likewise
humanFC_outputFile = paste0(humanFC_inputFile,".overlaps")
message("Intersect BED HumanFC:",appendLF=F)
system2('intersectBed', 
        paste('-wao -f 0.05 -a', paste0(dataDir,bin_outputFile,".bed"), '-b', paste0(dataDir,humanFC_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,humanFC_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

humanFC_inputFile = humanFC_outputFile
humanFC_outputFile = paste0(humanFC_inputFile,".dropped")
message("Drop columns in HumanFC:",appendLF=F)
system2('cut', 
        paste('-f1-4,9-297',paste0(dataDir,humanFC_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,humanFC_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

message("For multiple overlap with the same bin, pick the max-overlapped one: in HumanFC:",appendLF=F)
con <- file(paste0(dataDir,"mergedPeakHeightMatrix_HumanFC_filtered.bed"),"r")
header <- readLines(con,n=1) %>% strsplit("\t") %>% do.call(c,.)
close(con)
humanFC_inputFile = humanFC_outputFile
humanFC_outputFile = paste0(humanFC_inputFile,".filtered")
dat <- fread(paste0(dataDir,humanFC_inputFile,".bed"), sep="\t", header=F)
dat <- dat %>% group_by(V4) %>% slice(which.max(V293)) %>% select(-c(V293))
colnames(dat) <- header
fwrite(dat, file=paste0(dataDir,humanFC_outputFile,".bed"), sep="\t")
message("Done",appendLF=T)

message("Replace dots with 0 (comes from intersectBed -wao): in HumanFC:",appendLF=F)
humanFC_inputFile = humanFC_outputFile
humanFC_outputFile = paste0(humanFC_inputFile,".fixed")
system2('sed', 
        paste('\'s/\\./0/g\'',paste0(dataDir,humanFC_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,humanFC_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

message("Sort bins: in HumanFC:",appendLF=F)
humanFC_inputFile = humanFC_outputFile
humanFC_outputFile = paste0(humanFC_inputFile,".sorted")
system2('sort', 
        paste('-k 1,1 -k2,2n',paste0(dataDir,humanFC_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,humanFC_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

# ---------------------------------------------- EpiMap
EpiMap_inputFile = "mergedPeakHeightMatrix_EpiMap_filtered"
EpiMap_outputFile = paste0(EpiMap_inputFile,".overlaps")
message("Intersect BED EpiMap:",appendLF=F)
system2('intersectBed', 
        paste('-wao -f 0.05 -a', paste0(dataDir,bin_outputFile,".bed"), '-b', paste0(dataDir,EpiMap_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,EpiMap_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

EpiMap_inputFile = EpiMap_outputFile
EpiMap_outputFile = paste0(EpiMap_inputFile,".dropped")
message("Drop columns in EpiMap:",appendLF=F)
system2('cut', 
        paste('-f1-4,9-159',paste0(dataDir,EpiMap_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,EpiMap_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

message("For multiple overlap with the same bin, pick the max-overlapped one: in EpiMap:",appendLF=F)
con <- file(paste0(dataDir,"mergedPeakHeightMatrix_EpiMap_filtered.bed"),"r")
header <- readLines(con,n=1) %>% strsplit("\t") %>% do.call(c,.)
close(con)
EpiMap_inputFile = EpiMap_outputFile
EpiMap_outputFile = paste0(EpiMap_inputFile,".filtered")
dat <- fread(paste0(dataDir,EpiMap_inputFile,".bed"), sep="\t", header=F)
dat <- dat %>% group_by(V4) %>% slice(which.max(V155)) %>% select(-c(V155))
colnames(dat) <- header
fwrite(dat, file=paste0(dataDir,EpiMap_outputFile,".bed"), sep="\t")
message("Done",appendLF=T)

message("Replace dots with 0 (comes from intersectBed -wao): in HumanFC:",appendLF=F)
EpiMap_inputFile = EpiMap_outputFile
EpiMap_outputFile = paste0(EpiMap_inputFile,".fixed")
system2('sed', 
        paste('\'s/\\./0/g\'',paste0(dataDir,EpiMap_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,EpiMap_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

message("Sort bins: in EpiMap:",appendLF=F)
EpiMap_inputFile = EpiMap_outputFile
EpiMap_outputFile = paste0(EpiMap_inputFile,".sorted")
system2('sort', 
        paste('-k 1,1 -k2,2n',paste0(dataDir,EpiMap_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,EpiMap_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

# ---------------------------------------------- ENCODE_DCC TFs
TFs_inputFile = "final.tf"
TFs_outputFile = paste0(TFs_inputFile,".overlaps")
message("Intersect BED TFs:",appendLF=F)
system2('intersectBed', 
        paste('-wao -f 0.05 -a', paste0(dataDir,bin_outputFile,".bed"), '-b', paste0(dataDir,TFs_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,TFs_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

TFs_inputFile = TFs_outputFile
TFs_outputFile = paste0(TFs_inputFile,".dropped")
message("Drop columns in TFs:",appendLF=F)
system2('cut', 
        paste('-f1-4,9-137',paste0(dataDir,TFs_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,TFs_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

message("For multiple overlap with the same bin, pick the max-overlapped one: in TFs:",appendLF=F)
con <- file(paste0(dataDir,"final.tf.bed"),"r")
header <- readLines(con,n=1) %>% strsplit("\t") %>% do.call(c,.)
close(con)
TFs_inputFile = TFs_outputFile
TFs_outputFile = paste0(TFs_inputFile,".filtered")
dat <- fread(paste0(dataDir,TFs_inputFile,".bed"), sep="\t", header=F)
dat <- dat %>% group_by(V4) %>% slice(which.max(V133)) %>% select(-c(V133))
colnames(dat) <- header
fwrite(dat, file=paste0(dataDir,TFs_outputFile,".bed"), sep="\t")
message("Done",appendLF=T)

message("Replace dots with 0 (comes from intersectBed -wao): in TFs:",appendLF=F)
TFs_inputFile = TFs_outputFile
TFs_outputFile = paste0(TFs_inputFile,".fixed")
system2('sed', 
        paste('\'s/\\./0/g\'',paste0(dataDir,TFs_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,TFs_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)

message("Sort bins: in TFs:",appendLF=F)
TFs_inputFile = TFs_outputFile
TFs_outputFile = paste0(TFs_inputFile,".sorted")
system2('sort', 
        paste('-k 1,1 -k2,2n',paste0(dataDir,TFs_inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,TFs_outputFile,".bed"), 
        wait=T)
message("Done",appendLF=T)


# non-Zero bins extraction in Shell (using awk): three file outputs
message("Getting Non-zero genomic bins:",appendLF=F)
# HumanFC
inputFile = "mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.filtered.fixed.sorted"
outputFile = "HumanFC_nonZero.binInfo"
system2('awk', 
        paste('-F \'\t\' \' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }\'',paste0(dataDir,inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,outputFile,".bed"), 
        wait=T)
# EpiMap
inputFile = "mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.filtered.fixed.sorted"
outputFile = "EpiMap_nonZero.binInfo"
system2('awk', 
        paste('-F \'\t\' \' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }\'',paste0(dataDir,inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,outputFile,".bed"), 
        wait=T)
# TFs
inputFile = "final.tf.overlaps.dropped.filtered.fixed.sorted"
outputFile = "ENCODE_nonZero.binInfo"
system2('awk', 
        paste('-F \'\t\' \' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }\'',paste0(dataDir,inputFile,".bed"), sep=' '), 
        stdout=paste0(dataDir,outputFile,".bed"), 
        wait=T)

message("Done",appendLF=T)

# Union of non-Zero bins: one file output
setwd(dataDir)
unionFileOutput = "HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed"
# binIDs got damaged somehow (ie. scientific notation appears) - don't know when and why, so need to reconstruct
epi <- read.table("EpiMap_nonZero.binInfo.bed", sep='\t', header=F); epi <- cbind(epi[,-4],paste0(epi[,1],"_",epi[,2],"_",epi[,3]))
human <- read.table("HumanFC_nonZero.binInfo.bed", sep='\t', header=F);  human <- cbind(human[,-4],paste0(human[,1],"_",human[,2],"_",human[,3]))
tf <- read.table("ENCODE_nonZero.binInfo.bed", sep='\t', header=F); tf <- cbind(tf[,-4],paste0(tf[,1],"_",tf[,2],"_",tf[,3]))
colnames(human)=colnames(epi)=colnames(tf) <- c("chr","start","end","id")

# perform Union of records (bininfo); Ignore the warnings (auto-coercing of columns is helpful here)
human.epi <- dplyr::union(human,epi)
human.epi.tf <- dplyr::union(human.epi,tf)

fwrite(human.epi.tf,file=paste0(dataDir, unionFileOutput), sep="\t", row.names=F, quote=F)

# 

# Extract genomic bins and lables and combine in a single file: one File output

library(data.table)
library(dplyr)

library("BSgenome.Hsapiens.UCSC.hg19")
hg <- BSgenome.Hsapiens.UCSC.hg19
flankingLength <- flanking
inputFile <- unionFileOutput
DNAoutputFile <- "HumanFC_ENCODE_EpiMap_nonZero.bin.Seq"
nonZerobins <- fread(paste0(dataDir,inputFile), sep="\t", header=T)
seq <- getSeq(hg,nonZerobins$chr,start=nonZerobins$start-flankingLength,end=nonZerobins$end+flankingLength)
# seq <- tryCatch({
#   getSeq(hg,nonZerobins$chr,start=nonZerobins$start-flankingLength,end=nonZerobins$end+flankingLength)
# },error=function(e){
#   DNAString(paste(rep('N',2*flankingLength+binSize),collapse = ''))
# }
# )
seq <- as.character(as.data.frame(seq)[[1]])
nonZerobins.seq <- cbind(nonZerobins,seq)
colnames(nonZerobins.seq) <- c(colnames(nonZerobins),"dna.seq")
fwrite(nonZerobins.seq, file=paste0(dataDir, DNAoutputFile,".bed"), sep="\t", row.names=F, quote=F)

message("Extract labels for non-zero bins from each data files (HumanFC, EpiMap and ENCODE_TFs):",appendLF=F)
# HUMANFC
inputFile1 = "HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed"
inputFile2 = "mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.filtered.fixed.sorted.bed"
HumanFC_nonzero_outputFile = "HumanFC_nonzero_labels"
system2('awk', 
        paste0('-F \"\t\" \'FILENAME==\"',
               paste0(dataDir,inputFile1),
               '\"{A[$1$2$3]=$1$2$3} FILENAME==\"',
               paste0(dataDir,inputFile2),
               '\"{if(A[$1$2$3]==$1$2$3){print}}\' ',
               paste0(dataDir,inputFile1),' ',
               paste0(dataDir,inputFile2)), 
        stdout=paste0(dataDir,HumanFC_nonzero_outputFile,".bed"), 
        wait=T)

# EpiMap
inputFile1 = "HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed"
inputFile2 = "mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.filtered.fixed.sorted.bed"
EpiMap_nonzero_outputFile = "EpiMap_nonzero_labels"
system2('awk', 
        paste0('-F \"\t\" \'FILENAME==\"',
               paste0(dataDir,inputFile1),
               '\"{A[$1$2$3]=$1$2$3} FILENAME==\"',
               paste0(dataDir,inputFile2),
               '\"{if(A[$1$2$3]==$1$2$3){print}}\' ',
               paste0(dataDir,inputFile1),' ',
               paste0(dataDir,inputFile2)), 
        stdout=paste0(dataDir,EpiMap_nonzero_outputFile,".bed"), 
        wait=T)

# TFs
inputFile1 = "HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed"
inputFile2 = "final.tf.overlaps.dropped.filtered.fixed.sorted.bed"
TFs_nonzero_outputFile = "ENCODE_TFs_nonzero_labels"
system2('awk', 
        paste0('-F \"\t\" \'FILENAME==\"',
               paste0(dataDir,inputFile1),
               '\"{A[$1$2$3]=$1$2$3} FILENAME==\"',
               paste0(dataDir,inputFile2),
               '\"{if(A[$1$2$3]==$1$2$3){print}}\' ',
               paste0(dataDir,inputFile1),' ',
               paste0(dataDir,inputFile2)), 
        stdout=paste0(dataDir,TFs_nonzero_outputFile,".bed"), 
        wait=T)

message("Done",appendLF=T)

# HumanFC_nonzero_outputFile = "HumanFC_nonzero_labels"
# EpiMap_nonzero_outputFile = "EpiMap_nonzero_labels"
# TFs_nonzero_outputFile = "ENCODE_TFs_nonzero_labels"
message("Print all the nonzero bins with sequence and labels:",appendLF=F)
outputFile <- "HumanFC_ENCODE_EpiMap_nonZero.bin.Seq_Labels"
dna.dat <- fread(paste0(dataDir, DNAoutputFile, ".bed"), sep="\t", header=T)
human <- fread(paste0(dataDir, HumanFC_nonzero_outputFile, ".bed"), sep="\t", header=T)
human$id <- paste0(human$chr, "_", human$start, "_", human$end)
epi <- fread(paste0(dataDir, EpiMap_nonzero_outputFile, ".bed"), sep="\t", header=T)
epi$id <- paste0(epi$chr, "_", epi$start, "_", epi$end)
tf <- fread(paste0(dataDir, TFs_nonzero_outputFile, ".bed"), sep="\t", header=T)
tf$id <- paste0(tf$chr, "_", tf$start, "_", tf$end)

output <- cbind(human, epi[which(epi$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
output_full <- cbind(dna.dat[which(dna.dat$id %in% human$id), ], human[,-c(1:4)], epi[which(epi$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
colnames(output) <- c(colnames(human), colnames(epi)[-c(1:4)], colnames(tf)[-c(1:4)])
colnames(output_full) <- c(colnames(dna.dat), colnames(human)[-c(1:4)], colnames(epi)[-c(1:4)], colnames(tf)[-c(1:4)])
fwrite(output,file="HumanFC_ENCODE_EpiMap_nonZero.bin.Labels.bed", sep="\t", row.names=F, quote=F)
fwrite(output_full,file=paste0(dataDir, outputFile, ".bed"), sep="\t", row.names=F, quote=F)

message("Done",appendLF=T)

# split (train and test (default: 'chr')) and numpy in Shell (python code): four File outputs
# remember to explicitely add run permission to the python script (chmod +x split_and_Numpy_V2.py)
system2("python", 
        paste(paste0(dataDir,"split_and_Numpy_V2.py"), 
              "--datadir", dataDir,
              "--datafilename HumanFC_ENCODE_EpiMap_nonZero.bin.Seq_Labels",
              "--valid_chr_id 1", 
              sep = ' '), 
        wait=T)

# Train and Test with DL model (5 epochs) in Shell and output a single accuracy (median) of all labels from final epoch's vallidation
# Use Google CoLab

