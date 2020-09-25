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
