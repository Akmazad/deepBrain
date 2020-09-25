## Peak filtering
We selected rows that have at least 50% samples have a signal "1" (defined as peak height > 0) using [```peak_filtering_single_Sample.R```](https://github.com/Akmazad/deepPsych/blob/master/Data%20Processing/Psychencode_June2019/Single%20Label-based%20exp/peak_filtering_single_Sample.R)) script. It outputs filtered data in BED format along with binarized peak value at threshold 0.

|Name|nSample (only Controls)|nPeaks (total)| nPeaks (>50% samples have non-zero peak-height) |
|---|---|---|---|
|HumanFC|137 (out of 288)|197,263|15,867|

## Bin overlapping
Next, we need to intersect each peaks with chromosomal bins of fixed-width. Following commands were used to get it done.

```sh
intersectBed -wao -f 0.05 -a hg19_bins_200bp.bed -b mergedPeakHeightMatrix_HumanFC_filtered_single_label.bed > mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.bed
```

## Post-processing
- We need to drop few information that aren't relevant (comes from peaks' binIDs after [```intersectBed -wao```](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)).
```sh
cut -f1-4,9-10 mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.bed > mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.dropped.bed
```

- For the same bin that overlaps with multiple peak vectors, we should chose the one with max overlap, i.e. the last column indicates overlap ammount after running [```intersectBed -wao```](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)).

```r
library(dplyr)
library(data.table)
setwd("/srv/scratch/z3526914/DeepBrain/Data")
# for HumanFC
# read the header (i.e. sample names)
con <- file("mergedPeakHeightMatrix_HumanFC_filtered_single_label.bed","r")
header <- readLines(con,n=1) %>% strsplit("\t") %>% do.call(c,.)
close(con)
dat <- fread("mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.dropped.bed", sep="\t", header=F)
dat <- dat %>% group_by(V4) %>% slice(which.max(V6)) %>% select(-c(V6))
colnames(dat) <- header
fwrite(dat, file="mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.dropped.filtered.bed", sep="\t")

```

- Replace all the dots (comes from the [```intersectBed -wao```](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)) when no matches are found.
```sh
sed 's/\./0/g' mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.dropped.filtered.bed > mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.dropped.fixed.filtered.bed
```

## Sorting bins
This subsection sorts bins (they are in BED format) by chromosome then by start position (same as [```this subsection```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/README.md#28-sorting-bins)).
```sh
sort -k 1,1 -k2,2n mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.dropped.fixed.filtered.bed > mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.dropped.fixed.filtered.sorted.bed
```

## 3. (Non-zero genomic bins) based pipeline
- Genomic Bins (sized = 200bp) with at least one signal (binary 1) found among all samples
- Scripts used:
```sh
# for saving non-zero binInfo
awk -F '\t' ' {for(i=5; i<=NF; i++) if ($i != -1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5; break;} }' mergedPeakHeightMatrix_HumanFC_filtered_single_label.overlaps.dropped.fixed.filtered.sorted.bed > HumanFC_single_label_nonZero.binInfo.bed
```
## fix the bins ending when it goes beyond the refWidth

## Extract DNA sequences
```r
setwd("/srv/scratch/z3526914/DeepBrain/Data")
library(data.table)
library(dplyr)

library("BSgenome.Hsapiens.UCSC.hg19")
hg <- BSgenome.Hsapiens.UCSC.hg19

flankingLength = 500
bins <- fread("HumanFC_single_label.binInfo.bed", sep="\t", header=T)
# fix any bin ending+flankingWidth is beyond the length of the chromosome
chr.refWidth = fread("hg19.chrom.sizes.txt") %>% dplyr::rename(chr = V1, refWidth = V2)
bins <- dplyr::left_join(bins, chr.refWidth, by = "chr")
bins <- bins %>% dplyr::filter(end + flankingLength < refWidth)  %>% dplyr::select(-refWidth)

seq <- getSeq(hg, bins$chr, start = bins$start - flankingLength, end = bins$end + flankingLength)
seq <- as.character(as.data.frame(seq)[[1]])
bins.seq <- cbind(bins,seq)
colnames(bins.seq) <- c(colnames(bins),"dna.seq")
fwrite(bins.seq, file="HumanFC_single_label.bin.Seq.bed", sep="\t", row.names=F, quote=F)
```

################## End of data-processing Pipeline ##############


