# Introduction
In this data collection, we've gathered and preprocessed (using [```peak_processing_new.R```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/peak_processing_new.R)) CHIPseq and ATACseq data from Psychencode. Details to follow.
# Data
 - EpiMap: 150 samples; 479,476 peaks
 - HumanFC: 288 samples; 197,263 peaks

# Method
## Peak filtering
We dropped rows for which there are no more than 1 sample peaks found by using [```peak_filtering.R```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/Psychencode_June2019/peak_filtering.R)) script.
# Result
