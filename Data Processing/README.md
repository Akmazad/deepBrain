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
|HumanFC|490,23|118,347|
|ENCODE TFs|2,441,723|725,276|



