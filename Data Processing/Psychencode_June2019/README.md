# Introduction
In this data collection, we've gathered and preprocessed (using [```peak_processing_new.R```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/peak_processing_new.R)) CHIPseq and ATACseq data from Psychencode. Details to follow.
# Data
|Name|nSample|nPeaks|
|---|---|---|
|EpiMap|150|479,476|
|HumanFC|288 samples|197,263 peaks|

# Method
## Peak filtering
We dropped rows for which there are no more than 1 sample peaks found by using [```peak_filtering.R```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/peak_filtering.R)) script. It outputs filtered data in BED format along with binarized peak value yielding 353.
# Result
