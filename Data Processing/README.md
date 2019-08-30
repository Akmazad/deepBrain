# Data preprocessing for Deep learning training by analysing RNA-seq, CHIP-seq, ATAC-seq, and Transcription factor peaks (also CHIP-seq from UCSC-DCC)
Our pipeline considers only those chromosomal bins for DL training that has at least one TF features activated. Hence, we first constructed bins with at least one TF signal found, and after that augmented other features from e.g. EpiMap, HumanFC, or other sources.

## Processed data
### Merged and filtered peaks
- Filtering thresholds applied: peak value > 0; nSample >= 2

|Name|nSample|nPeaks|
|---|---|---|
|EpiMap|150|353,566|
|HumanFC|288|118,347|
|ENCODE TFs|128|725,276|

### Non-zero genomic bins
- Genomic Bins (sized = 200bp) with at least one signal (binary 1) found among all samples

|Name|nBins|nPeak coordinates (filtered)|
|---|---|---|
|EpiMap|1,744,883|353,566|
|HumanFC|490,233|118,347|
|ENCODE TFs|2,441,723|725,276|
|union|3,528,533|---|
- Scripts used:
```sh
# for saving non-zero binInfo
awk -F '\t' ' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }' mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.fixed.filtered.sorted.dat > EpiMap_nonZero.binInfo.dat
awk -F '\t' ' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }' mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.fixed.filtered.sorted.dat > HumanFC_nonZero.binInfo.dat
awk -F '\t' ' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }' final.dat.tf.overlaps.dropped.fixed.filtered.sorted.dat > ENCODE_nonZero.binInfo.dat
```
```r
# for making Union of all non-zero binInfo
# on KATANA (head node): 
setwd('/srv/scratch/z3526914/DeepBrain/Data')
library(dplyr)
library(data.table)

# binIDs got damaged somehow (ie. scientific notation appears) - don't know when and why, so need to reconstruct
epi <- read.table("EpiMap_nonZero.binInfo.dat", sep='\t', header=F); epi <- cbind(epi[,-4],paste0(epi[,1],"_",epi[,2],"_",epi[,3]))
human <- read.table("HumanFC_nonZero.binInfo.dat", sep='\t', header=F);  human <- cbind(human[,-4],paste0(human[,1],"_",human[,2],"_",human[,3]))
tf <- read.table("ENCODE_nonZero.binInfo.dat", sep='\t', header=F); tf <- cbind(tf[,-4],paste0(tf[,1],"_",tf[,2],"_",tf[,3]))
colnames(human)=colnames(epi)=colnames(tf) <- c("chr","start","end","id")

# perform Union of records (bininfo); Ignore the warnings (auto-coercing of columns is helpful here)
human.epi <- dplyr::union(human,epi)
human.epi.tf <- dplyr::union(human.epi,tf)
# nrow(human.epi.tf):
# [1] 3528533
fwrite(human.epi.tf,file="HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.dat", sep="\t", row.names=F, quote=F)
```

## Extract genomic data (dna seq) for non-zero bins
Need to run on KATANA ([```ExtractDNAseq_KATANA.sh```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/ExtractDNAseq_KATANA.sh))
```r
rm(list = ls(all.names = TRUE))
setwd('/srv/scratch/z3526914/DeepBrain/Data')
library(data.table)
library(dplyr)
# read binned DNA seq (human)
hg <- fread("hg19.binwise.fasta.200bp.bed", sep="\t", header=F)
colnames(hg) <- c("chr","start","end","dna.seq","id","strand")

# read non-zero bins
nonZerobins <- fread("HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.dat", sep="\t", header=T)
nonZerobins.seq <- cbind(nonZerobins,hg[which(hg$id %in% nonZerobins$id), "dna.seq"])
colnames(nonZerobins.seq) <- c(colnames(nonZerobins),"dna.seq")
fwrite(nonZerobins.seq, file="HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.dat", sep="\t", row.names=F, quote=F)
```

## Extract Labels (binary signals) for non-zero bins

