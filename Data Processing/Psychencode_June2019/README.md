# Introduction
In this data collection, we've gathered and preprocessed (using [```peak_processing_new.R```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/peak_processing_new.R)) CHIPseq and ATACseq data from Psychencode. Details to follow.
# Data
|Name|nSample|nPeaks|
|---|---|---|
|EpiMap|150|479,476|
|HumanFC|288|197,263|

# Method
## Peak filtering
We dropped rows for which there are no more than 1 sample peaks found by using [```peak_filtering.R```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/peak_filtering.R)) script. It outputs filtered data in BED format along with binarized peak value at threshold 0.

|Name|nSample|nPeaks|
|---|---|---|
|EpiMap|150|353,566|
|HumanFC|288|118,347|

This filtered peak matrices can be augmented with TF features sets using [```Merge_binned_TF_profiles_with_other_features.sh```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/RNA-seq/Merge_binned_TF_profiles_with_other_features.sh) script on Raijin *hugemem* queue.

## Bin overlapping
Next, we need to intersect each peaks with chromosomal bins of fixed-width. Following commands were used to get it done.

```shell
intersectBed -wao -f 0.05 -a hg19_bins_200bp.bed -b mergedPeakHeightMatrix_HumanFC_filtered.bed > mergedPeakHeightMatr.overlaps.bed
intersectBed -wao -f 0.05 -a hg19_bins_200bp.bed -b mergedPeakHeightMatrix_EpiMap_filtered.bed > mergedPeakHeightMatrix_EpiMap_filtered.overlaps.bed
```
We need to drop few information that aren't relevant (comes from peaks' binIDs after ```intersectBed```).
```awk
awk '{$5=$6=$7=$8=$9=$298=""; print $0}' mergedPeakHeightMatrix_HumanFC_filtered.overlaps.bed > mergedPeakHeightMatrix_HumanFC_filtered.overlaps.dropped.bed
awk '{$5=$6=$7=$8=$9=$298=""; print $0}' mergedPeakHeightMatrix_EpiMap_filtered.overlaps.bed > mergedPeakHeightMatrix_EpiMap_filtered.overlaps.dropped.bed
```
# Result
