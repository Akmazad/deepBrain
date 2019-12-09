# Introduction
In this data collection, we've gathered and preprocessed (using [```peak_processing_new.R```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/peak_processing_new.R)) CHIPseq and ATACseq data from Psychencode. Details to follow.
# Data
|Name|nSample|nPeaks|
|---|---|---|
|EpiMap|150|479,476|
|HumanFC|288|197,263|
|CAGE|-(Jerry to fill)-|-(Jerry to fill)-|

# Method
## Data Collection [Contributor: Jerry Offor]
EpiMap and HumanFC studies were download from the PsychEncode Consortium using synapseClient(https://www.synapse.org). All samples from EpiMap came from NIMH Human Brain Collection core (HBCC).After dissection, they were shipped for sample preparation and CHIPseq to Ichan School of Medicince-Mt Sinai (ISMMS). Samples from HumanFc study came from Mount Sinai NIH Brain Bank and Tissue Repository and then shipped to Duke for preparation of samples and ATAC-Seq .

### EpiMap Data-ChiPseq
BWA (VERSION 0.7.8) was used to align paired end FASTQ files to the human genome (HG19). Duplicates were mark into bamfiles using Picard (version 1.1.12) MarkDups. Filtering of improperly paired reads and multi-mapped reads was done with Samtools with defined parameters. The peaks (Narrow peaks) from Duplicate marked were call with Macs2.

### HumanFC Data-ATAC_seq
Fastq files were first processed with cut-adaptor (version 1.2.0) to remove low quality read and adaptors and then aligned with bowtie 2(version 2.1.0) to the human genome (hg19) with default parameters. All reads during alignment were treated as single read sequence even if sequencing libraries were single or paired end. Samtools (version 0.1.18) was used to sort bam files, picard.jar MarkDuplicate was used to filter-out duplicates and then converted to bed format using bamToBed before narrow peaks were called using MACS2 with a threshold FDR < 0.01.

Then Bedtools were used to combine all bed files

The table below shows the number of samples per study selected from PsychEncode Consortium:
[TBD: from Jerry's Email]

## Peak filtering
We dropped rows for which there are no more than 1 sample peaks found by using [```peak_filtering.R```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/peak_filtering.R)) script. It outputs filtered data in BED format along with binarized peak value at threshold 0.

|Name|nSample|nPeaks|
|---|---|---|
|EpiMap|150|353,566|
|HumanFC|288|118,347|
|CAGE|136|163,018|

This filtered peak matrices can be augmented with TF features sets using [```Merge_binned_TF_profiles_with_other_features.sh```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/RNA-seq/Merge_binned_TF_profiles_with_other_features.sh) script on Raijin *hugemem* queue.

## Bin overlapping
Next, we need to intersect each peaks with chromosomal bins of fixed-width. Following commands were used to get it done.

```sh
intersectBed -wao -f 0.05 -a hg19_bins_200bp.bed -b mergedPeakHeightMatrix_HumanFC_filtered.bed > mergedPeakHeightMatrix_HumanFC_filtered.overlaps.bed
intersectBed -wao -f 0.05 -a hg19_bins_200bp.bed -b mergedPeakHeightMatrix_EpiMap_filtered.bed > mergedPeakHeightMatrix_EpiMap_filtered.overlaps.bed
intersectBed -wao -f 0.05 -a hg19_bins_200bp.bed -b Brain_CagePeaks_filtered.BED > Brain_CagePeaks_filtered.overlaps.bed
```
## Post-processing
- We need to drop few information that aren't relevant (comes from peaks' binIDs after [```intersectBed -wao```](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)).
```sh
cut -f1-4,10-298 mergedPeakHeightMatrix_HumanFC_filtered.overlaps.bed > mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.bed
cut -f1-4,10-160 mergedPeakHeightMatrix_EpiMap_filtered.overlaps.bed > mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.bed
cut -f1-4,10-146 Brain_CagePeaks_filtered.overlaps.bed > Brain_CagePeaks_filtered.overlaps.dropped.bed
```
- For the same bin that overlaps with multiple peak vectors, we should chose the one with max overlap, i.e. the last column indicates overlap ammount after running [```intersectBed -wao```](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)). This subsection follows similar steps in [```TF processing pipeline```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/README.md#27-filter-similar-overlapping-bins-with-the-max-overlap-size-last-column). But the codes are copied here though. NOTE: THIS SCRIPT FOR HumanFC EXHAUSTS RNA MACHINE'S MEMORY: HENCE, KATANA IS APPLIED ([```HumanFC_post_processing.sh```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/HumanFC_post_processing.sh)).
```r
library(dplyr)
library(data.table)
# for HumanFC
# read the header (i.e. sample names)
con <- file("mergedPeakHeightMatrix_HumanFC_filtered.bed","r")
header <- readLines(con,n=1) %>% strsplit("\t") %>% do.call(c,.)
close(con)
dat <- fread("mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.bed", sep="\t", header=F)
dat <- dat %>% group_by(V4) %>% slice(which.max(V293)) %>% select(-c(V293))
colnames(dat) <- header
fwrite(dat, file="mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.filtered.bed", sep="\t")

# for EpiMap
# read the header (i.e. sample names)
con <- file("mergedPeakHeightMatrix_EpiMap_filtered.bed","r")
header <- readLines(con,n=1) %>% strsplit("\t") %>% do.call(c,.)
close(con)
dat <- fread("mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.bed", sep="\t", header=F)
dat <- dat %>% group_by(V4) %>% slice(which.max(V155)) %>% select(-c(V155))
colnames(dat) <- header
fwrite(dat, file="mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.filtered.bed", sep="\t")

# for CAGE
# read the header (i.e. sample names)
con <- file("Brain_CagePeaks_filtered.bed","r")
header <- readLines(con,n=1) %>% strsplit("\t") %>% do.call(c,.)
close(con)
dat <- fread("Brain_CagePeaks_filtered.overlaps.dropped.bed", sep="\t", header=F)
dat <- dat %>% group_by(V4) %>% slice(which.max(V141)) %>% select(-c(V141))
colnames(dat) <- header
fwrite(dat, file="Brain_CagePeaks_filtered.overlaps.dropped.filtered.bed", sep="\t")

```

- Replace all the dots (comes from the [```intersectBed -wao```](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)) when no matches are found.
```sh
sed 's/\./0/g' mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.filtered.bed > mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.fixed.filtered.bed
sed 's/\./0/g' mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.filtered.bed > mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.fixed.filtered.bed
sed 's/\./0/g' Brain_CagePeaks_filtered.overlaps.dropped.filtered.bed > Brain_CagePeaks_filtered.overlaps.dropped.filtered.fixed.filtered.bed
```

## Sorting bins
This subsection sorts bins (they are in BED format) by chromosome then by start position (same as [```this subsection```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/README.md#28-sorting-bins)).
```sh
sort -k 1,1 -k2,2n mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.fixed.filtered.bed > mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.fixed.filtered.sorted.bed
sort -k 1,1 -k2,2n mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.fixed.filtered.bed > mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.fixed.filtered.sorted.bed
sort -k 1,1 -k2,2n Brain_CagePeaks_filtered.overlaps.dropped.filtered.fixed.filtered.bed > Brain_CagePeaks_filtered.overlaps.dropped.filtered.fixed.filtered.sorted.bed
```
################## End of data-processing Pipeline ##############

# Result




1) Few commands used for intermmediate sanity checks:
- counts number of unique binIDs in file:
```sh
cut -f 4 mergedPeakHeightMatrix_HumanFC_filtered.overlaps.bed | sort | uniq | wc -l
```
- print lines in a file where a string 'e+' (scientific notation) appears in the 4th column of a file:
```sh
awk '$4 ~ /e+/ { print }' final.dat.tf.overlaps.dropped.fixed.filtered.dat > freq.tf.bed
```
- count the number of times a word 'e+' apprears in the whole file
```sh
grep -o 'e+' mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.filtered.fixed.dat | wc -l
```
