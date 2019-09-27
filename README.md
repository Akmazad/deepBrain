# deepBrain
Deep learning in predicting functional effects of noncoding variants in brain-specific disorders

# Steps:
1. [```Data preprocessing```](https://github.com/Akmazad/deepBrain/tree/master/Data%20Processing) (e.g. RNA-seq, CHIP-seq, ATAC-seq, etc.)
2. [```Data```](https://github.com/Akmazad/deepBrain/tree/master/Data%20Processing#extract-genomic-data-dna-seq-for-non-zero-bins) and [```Label```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/README.md#extract-labels-binary-signals-for-non-zero-bins) extraction for genomic bins (of interest)
3. [```Numpy conversion```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/split_and_Numpy_V2.py) of data and labels (for DL usage)
4. DL training/Validation/Testing (demonstrated in [```Training 3```](https://github.com/Akmazad/deepBrain/tree/master/Training%203%20%5Bwith%20New%20data%5D) subdirectory)

Note: A complete pipeline (R-based) automating the steps 1-3 can be found in [```full_pipeline.R```](https://github.com/Akmazad/deepBrain/blob/master/Full%20Pipeline/full_pipeline.R).





