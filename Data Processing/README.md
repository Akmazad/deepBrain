# Data preprocessing for Deep learning training by analysing RNA-seq, CHIP-seq, ATAC-seq, and Transcription factor peaks (also CHIP-seq from UCSC-DCC)
Our pipeline considers only those chromosomal bins for DL training that has at least feature found. 

## 1. Processing pipeline
|Data|Used for|Pipeline Documentations|
|---|---|---|
|EpiMap (Chip-seq)|Sample peak Label extraction|[```ReadMe```](https://github.com/Akmazad/deepBrain/tree/master/Data%20Processing/Psychencode_June2019/README.md)|
|HumanFC (ATAC-seq)|Sample peak Label extraction|[```ReadMe```](https://github.com/Akmazad/deepBrain/tree/master/Data%20Processing/Psychencode_June2019/README.md)|
|ENCODE TFs (RNA-seq + ENCODE DCC)|TF label extraction|[```ReadMe```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/RNA-seq/README.md)|


## 2. Processed data
### 2.1 Merged and filtered peaks
- Filtering thresholds applied: peak value > 0; nSample >= 2

|Name|nSample|nPeaks|FileName|FileLocation|FileSize|
|---|---|---|---|---|---|
|EpiMap|150|353,566|mergedPeakHeightMatrix_EpiMap_filtered.bed|/Volumes/Data1/PROJECTS/DeepLearning/Test|123,333,174 byte|
|HumanFC|288|118,347|mergedPeakHeightMatrix_HumanFC_filtered.bed|/Volumes/Data1/PROJECTS/DeepLearning/Test|73,943,685 byte|
|ENCODE TFs|128|725,276|final.tf.bed|/Volumes/Data1/PROJECTS/DeepLearning/Test|220,373,843 byte|

## 3. (Non-zero genomic bins) based pipeline
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
awk -F '\t' ' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }' mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.fixed.filtered.sorted.bed > EpiMap_nonZero.binInfo.bed
awk -F '\t' ' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }' mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.fixed.filtered.sorted.bed > HumanFC_nonZero.binInfo.bed
awk -F '\t' ' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }' final.tf.overlaps.dropped.fixed.filtered.sorted.bed > ENCODE_nonZero.binInfo.bed
```
```r
# for making Union of all non-zero binInfo
# on KATANA (head node): 
setwd('/srv/scratch/z3526914/DeepBrain/Data')
library(dplyr)
library(data.table)

# binIDs got damaged somehow (ie. scientific notation appears) - don't know when and why, so need to reconstruct
epi <- read.table("EpiMap_nonZero.binInfo.bed", sep='\t', header=F); epi <- cbind(epi[,-4],paste0(epi[,1],"_",epi[,2],"_",epi[,3]))
human <- read.table("HumanFC_nonZero.binInfo.bed", sep='\t', header=F);  human <- cbind(human[,-4],paste0(human[,1],"_",human[,2],"_",human[,3]))
tf <- read.table("ENCODE_nonZero.binInfo.bed", sep='\t', header=F); tf <- cbind(tf[,-4],paste0(tf[,1],"_",tf[,2],"_",tf[,3]))
colnames(human)=colnames(epi)=colnames(tf) <- c("chr","start","end","id")

# perform Union of records (bininfo); Ignore the warnings (auto-coercing of columns is helpful here)
human.epi <- dplyr::union(human,epi)
human.epi.tf <- dplyr::union(human.epi,tf)
# nrow(human.epi.tf):
# [1] 3528533
fwrite(human.epi.tf,file="HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed", sep="\t", row.names=F, quote=F)
```

### 3.1 Extract genomic data (dna seq) for non-zero bins
Need to run on KATANA ([```ExtractDNAseq_KATANA.sh```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/ExtractDNAseq_KATANA.sh)) with following command (excerpt from the bash script):
```sh
Rscript /srv/scratch/z3526914/DeepBrain/Scripts/ExtractDNAseq_KATANA.R \
	/srv/scratch/z3526914/DeepBrain/Data \
	400 \
	HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed \
	HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.bed
```
### 3.2 Extract Labels (binary signals) for non-zero bins
- Extract labels for non-zero bins from each data files (HumanFC, EpiMap and ENCODE_TFs)
```sh
# its checking binInfo of both files (chr, start and end coordinates of bins)
awk -F "\t" 'FILENAME=="HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed"{A[$1$2$3]=$1$2$3} FILENAME=="mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.fixed.filtered.sorted.bed"{if(A[$1$2$3]==$1$2$3){print}}' HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.fixed.filtered.sorted.bed > HumanFC_nonzero_labels.bed

awk -F "\t" 'FILENAME=="HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed"{A[$1$2$3]=$1$2$3} FILENAME=="mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.fixed.filtered.sorted.bed"{if(A[$1$2$3]==$1$2$3){print}}' HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.fixed.filtered.sorted.bed > EpiMap_nonzero_labels.bed

awk -F "\t" 'FILENAME=="HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed"{A[$1$2$3]=$1$2$3} FILENAME=="final.tf.overlaps.dropped.fixed.filtered.sorted.bed"{if(A[$1$2$3]==$1$2$3){print}}' HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed final.tf.overlaps.dropped.fixed.filtered.sorted.bed > ENCODE_TFs_nonzero_labels.bed

```
- Merge all labels. Need to run on KATANA ([```ExtractLabels_KATANA.sh```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/ExtractLabels_KATANA.sh)). DNA sequences (Data) will be also augmented
```r
rm(list = ls(all.names = TRUE))
setwd('/srv/scratch/z3526914/DeepBrain/Data')
library(data.table)
library(dplyr)

# Label data files: 
# 1. HumanFC_nonzero_labels.bed
# 2. EpiMap_nonzero_labels.bed
# 3. ENCODE_TFs_nonzero_labels.bed
# 4. HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.bed

dna.dat <- fread("HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.bed", sep="\t", header=T) # check ids: "grep -o 'e+' HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.bed | wc -l" should return 0
dna.dat$id = paste0(dna.dat$chr, "_", dna.dat$start, "_", dna.dat$end)

human <- fread("HumanFC_nonzero_labels.bed", sep="\t", header=T)
human$id <- paste0(human$chr, "_", human$start, "_", human$end) # fix the 746 ids (scientific notation appread!!)
epi <- fread("EpiMap_nonzero_labels.bed", sep="\t", header=T)
epi$id <- paste0(epi$chr, "_", epi$start, "_", epi$end) # fix the 746 ids (scientific notation appread!!)
tf <- fread("ENCODE_TFs_nonzero_labels.bed", sep="\t", header=T)
tf$id <- paste0(tf$chr, "_", tf$start, "_", tf$end) # fix the 746 ids (scientific notation appread!!)

output <- cbind(human, epi[which(epi$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
output_full <- cbind(dna.dat[which(dna.dat$id %in% human$id), ], human[,-c(1:4)], epi[which(epi$id %in% human$id),-c(1:4)], tf[which(tf$id %in% human$id),-c(1:4)])
colnames(output) <- c(colnames(human), colnames(epi)[-c(1:4)], colnames(tf)[-c(1:4)])
colnames(output_full) <- c(colnames(dna.dat), colnames(human)[-c(1:4)], colnames(epi)[-c(1:4)], colnames(tf)[-c(1:4)])
fwrite(output,file="HumanFC_ENCODE_EpiMap_nonZero.bin.Labels.bed", sep="\t", row.names=F, quote=F)
fwrite(output_full,file="HumanFC_ENCODE_EpiMap_nonZero.bin.Seq_Labels.bed", sep="\t", row.names=F, quote=F)
```

### 3.3 Results
Final set of data (before entering DL pipeline) stats are as follows:

|Type|Filename|Location|nBins|nLabels|
|---|---|---|---|---|
|Genomic DNA|HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.bed|/Volumes/Data1/PROJECTS/DeepLearning/Test|3,528,533|---|
|Binary Labels|HumanFC_ENCODE_EpiMap_nonZero.bin.Labels.bed|/Volumes/Data1/PROJECTS/DeepLearning/Test|3,528,533|566|
|Data + Labels|HumanFC_ENCODE_EpiMap_nonZero.bin.Seq_Labels.bed|/Volumes/Data1/PROJECTS/DeepLearning/Test|3,528,533|566|

## 4. (TF-specific genomic bins) based pipeline
- Genomic Bins (sized = 200bp) with at least TF signal (binary 1) found among ENCODE TF genes

|Name|nBins|nPeak coordinates (filtered)|
|---|---|---|
|ENCODE TFs|2,441,723|725,276|

### 4.1 Get genomic bins with non-Zero TF signals
```sh
# for saving non-zero binInfo
awk -F '\t' ' {for(i=5; i<=NF; i++) if ($i == 1) {print $1"\t"$2"\t"$3"\t"$4; break;} }' final.tf.overlaps.dropped.fixed.filtered.sorted.bed > ENCODE_nonZero.binInfo.bed

# for adding the column headers
echo -e "chr\tstart\tend\tid\n$(cat ENCODE_nonZero.binInfo.bed)" > ENCODE_nonZero.binInfo.bed
```
### 4.2 Extract genomic data (dna seq) for non-zero bins
Need to run on KATANA ([```ExtractDNAseq_KATANA.sh```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/ExtractDNAseq_KATANA.sh)) with following command (excerpt from the bash script):
```sh
Rscript /srv/scratch/z3526914/DeepBrain/Scripts/ExtractDNAseq_KATANA.R \
	/srv/scratch/z3526914/DeepBrain/Data \
	400 \
	ENCODE_nonZero.binInfo.bed \
	HumanFC_ENCODE_EpiMap_tf_specific.bin.Seq.bed
```
Note, the input file (args[4]) is the 'ENCODE_nonZero.binInfo.bed', and only those bins (non-zero TFs) were considered for all other datasets. Hence the output file name 'HumanFC_ENCODE_EpiMap_tf_specific.bin.Seq.bed'.

### 4.2 Extract Labels (binary signals) for tf-specific bins
- Extract labels for tf-specific bins from each data files (HumanFC, EpiMap and ENCODE_TFs)
```sh
# its checking binInfo of both files (chr, start and end coordinates of bins)
awk -F "\t" 'FILENAME=="ENCODE_nonZero.binInfo.bed"{A[$1$2$3]=$1$2$3} FILENAME=="mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.fixed.filtered.sorted.bed"{if(A[$1$2$3]==$1$2$3){print}}' ENCODE_nonZero.binInfo.bed mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.fixed.filtered.sorted.bed > HumanFC_tf_specific_labels.bed

awk -F "\t" 'FILENAME=="ENCODE_nonZero.binInfo.bed"{A[$1$2$3]=$1$2$3} FILENAME=="mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.fixed.filtered.sorted.bed"{if(A[$1$2$3]==$1$2$3){print}}' ENCODE_nonZero.binInfo.bed mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.fixed.filtered.sorted.bed > EpiMap_tf_specific_labels.bed

awk -F "\t" 'FILENAME=="ENCODE_nonZero.binInfo.bed"{A[$1$2$3]=$1$2$3} FILENAME=="final.tf.overlaps.dropped.fixed.filtered.sorted.bed"{if(A[$1$2$3]==$1$2$3){print}}' ENCODE_nonZero.binInfo.bed final.tf.overlaps.dropped.fixed.filtered.sorted.bed > ENCODE_TFs_tf_specific_labels.bed

```
- Merge all labels. Need to run on KATANA ([```ExtractLabels_KATANA.sh```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/ExtractLabels_KATANA.sh)). DNA sequences (Data) will be also augmented. Command excerpt from the bash script:
```sh
Rscript /srv/scratch/z3526914/DeepBrain/Scripts/ExtractLabels_KATANA.R \
	/srv/scratch/z3526914/DeepBrain/Data \
	HumanFC_tf_specific_labels.bed \
	EpiMap_tf_specific_labels.bed \
	ENCODE_TFs_tf_specific_labels.bed \
	HumanFC_ENCODE_EpiMap_tf_specific.bin.Seq.bed \
	HumanFC_ENCODE_EpiMap_tf_specific.bin.Labels.bed \
	HumanFC_ENCODE_EpiMap_tf_specific.bin.Seq_Labels.bed
```
