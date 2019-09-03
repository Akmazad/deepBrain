# DeepBrain Training [version 3]
- This is a PyTorch implementation of the DeepBrain project. This project aims to predict the functional effects of non-coding variants from sequence data.
- In this version [V3], we've used Essays from HumanFC, EpiMap and TF profiles from ENCODE DCC to train deep learning model.

## Requirements
- python 3.6 or higher
- Pytorch 1.0.1
- Numpy
- Scipy
- sklearn

## Overview
| File | Description |
| --- | --- |
| [```split_and_Numpy_V2.py```](https://github.com/Akmazad/deepBrain/blob/master/Data%20Processing/split_and_Numpy_V2.py) | Prepares data for running DeepBrain models; train-validation split based on a chromosome |
| [```DL_input_Test_HET_Katana.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%203%20%5Bwith%20New%20data%5D/DL_input_Test_HET_Katana.py) | DeepBrain with static Convnet (e.g. [DeepSEA](https://github.com/FunctionLab/selene/blob/master/models/deepsea.py) and [DeeperDeepSEA](https://github.com/FunctionLab/selene/blob/master/tutorials/quickstart_training/deeperdeepsea.py)) models |


## Usage
### Data Preparation
The preprocessing steps yields data that contains both input DNA seqeunce and corresponding label for all chromosomes combined in a single file. Note, each record in this file has at least one signal present (i.e. they are nonZero rows). The file header looks like the following:

| chr | start | end | dna.seq | id | strand | 288 HumanFC features | 150 EpiMap feature | 128 TF features |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

DeepBrain models considers "dna.seq" as input data and the feature values (binary) as known label. Hence, we need to extract them from above file and convert them to numpy ndarrays. Moreover, we wanted to leave one chromosome data out while training, so that we can test the model with that. Therefore, we have to do train-validation split on both the value and label data based on a single chromosome. This process runs on Katana GPUs (for details, see the pbs script: [```DL_input_Test_HET_Katana.sh```](https://github.com/Akmazad/deepBrain/blob/master/Training%203%20%5Bwith%20New%20data%5D/DL_input_Test_HET_Katana.sh))

### Training and Validation
Initially, we are trying DeepSEA model-like architectures for training, which may be followed by the DARTs or other static ConvNet models with varied parameters settings. Training the DeeperDeepSEA-like model using [```DL_model_test.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/DL_model_test.py) and [```deepbrain2_dist.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/deepbrain2_dist.py), runs in single-GPU on Google CoLab ([```DL_input_test_HET_CoLab_notebook.ipynb```](https://github.com/Akmazad/deepBrain/blob/master/%20Training%203%20%5Bwith%20New%20data%5D/DL_input_test_HET_CoLab_notebook.ipynb)) and Raijin GPUs, respectively. 

To run on CPU, use following command:
```
 python3 DL_model_test_HET_Katana.py  \
--name deepbrainStaticConvnet \
--DataDir H:\ \
--TrainingDataFile temp_HET_trainingData_chr1_value.npy \
--TrainingLabelFile temp_HET_trainingData_chr1_label_all.npy \
--TestingDataFile temp_HET_validationData_chr1_value.npy \
--TestingLabelFile temp_HET_validationData_chr1_label_all.npy \
--nEpochs 15 \
--BATCH_SIZE 128
```

### Accuracy measures and Log reporting
For this version, we've checked two accuracy measurements. In each training/testing iteration (for a mini-batch), we report "true prediction ratio" for each Feature columns and take median accross three Feature categories (e.g. Accetylation, RNA-seq, and TFs).we also report AUC scores in a similar manner. In addition, progress logs are reported in a log file within the current data directory and project name (e.g. deepbrainStaticConvnet) subdirectory.


## Data Description

| Study | Assay Type | Nsamples | Number of peaks/Transcripts |
| --- | --- | --- | --- |
| HumanFC | ATAC-seq | 288 | 197,263 |
| EpiMap | ChIP-seq | 150 | 479,476 |
| UCSD-YALE (ba9) | ChIP-seq | 81 | 56,503 |
| UCSD-YALE (ba41) | ChIP-seq | 66 | 56,503 |
| UCSD-YALE (baVermis) | ChIP-seq | 62 | 38,069 |
| UCSD-YALE | RNA-seq | 209 | 392,033 |

## Thresholding Acetylation data
### Methodology
Choose peaks that has non-zero values 
