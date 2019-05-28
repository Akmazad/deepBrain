# DeepBrain
This is a PyTorch implementation of the DeepBrain project. This project aims to predict the functional effects of non-coding variants from sequence data.

## Requirements
- python 3.6 or higher
- Pytorch 1.0.1
- Numpy
- Scipy
- sklearn

## Overview
| File | Description |
| --- | --- |
| [```split_and_Numpy.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/split_and_Numpy.py) | Prepares data for running DeepBrain models; train-validation split based on a chromosome |
| [```DL_model_test.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/DL_model_test.py) | DeepBrain with static Convnet (e.g. [DeepSEA](https://github.com/FunctionLab/selene/blob/master/models/deepsea.py) and [DeeperDeepSEA](https://github.com/FunctionLab/selene/blob/master/tutorials/quickstart_training/deeperdeepsea.py)) models |
| [```DL_model_test_DARTS.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/DL_model_test_DARTS.py) | [DARTs](https://github.com/quark0/darts) (Differentiable Architecture Search) implementation of the ```DL_model_test.py``` * | 
| [```DL_model_test_DeepSEA_data.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/DL_model_test_DeepSEA_data.py) | ```DL_model_test.py``` with DeepSEA data |
| [```deepbrain2_dist.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/deepbrain2_dist.py) | Distributed (Multi-processs, Multi-GPU implementation of ```DL_model_test.py``` |
| [```deepbrain2_dist_DeepSEA_data.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/deepbrain2_dist_DeepSEA_data.py) | ```deepbrain2_dist.py``` with DeepSEA data |

## Usage
### Data Preparation
The preprocessing steps yields data that contains both input DNA seqeunce and corresponding label for all chromosomes combined in a single file. Note, each record in this file has at least one TF signal present. The file header looks like the following:

| chr | start | end | dna.seq | id | strand | 2 Accetylation features | 1 RNA-seq feature | 128 TF features |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

DeepBrain models considers "dna.seq" as input data and the feature values (binary) as known label. Hence, we need to extract them from above file and convert them to numpy ndarrays. Moreover, we wanted to leave one chromosome data out while training, so that we can test the model with that. Therefore, we have to do train-validation split on both the value and label data based on a single chromosome. This process runs on raijin (for details, see the pbs script: [```DL_input_TrainValid_Split_and_Numpy.sh```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/pbs%20scripts/DL_input_TrainValid_Split_and_Numpy.sh))

### Training and Validation
Initially, we are trying DeepSEA model-like architectures for training, which may be followed by the DARTs or other static ConvNet models by varying different parameters. Training the DeeperDeepSEA-like model using [```DL_model_test.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/DL_model_test.py) and [```deepbrain2_dist.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/deepbrain2_dist.py), runs in CPU and Raijin GPUs, respectively. 

To run on CPU, use following command:
```
 python3 DL_model_test.py  \
--name deepbrainStaticConvnet \
--DataDir H:\ \
--TrainingDataFile tempTrain5_trainingData_TF_filtered_chr5_value.npy \
--TrainingLabelFile tempTrain5_trainingData_TF_filtered_chr5_label_all.npy \
--TestingDataFile tempTrain5_validationData_TF_filtered_chr5_value.npy \
--TestingLabelFile tempTrain5_validationData_TF_filtered_chr5_label_all.npy \
--nEpochs 15 \
--BATCH_SIZE 128
```

For running [```deepbrain2_dist.py```](https://github.com/Akmazad/deepBrain/blob/master/Training%202/deepbrain2_dist.py) on Raijin GPUs on a distributed manner, please see the PBS script: [DeepBrain2_training_dist.sh](https://github.com/Akmazad/deepBrain/blob/master/Training%202/pbs%20scripts/DeepBrain2_training_dist.sh).

### Accuracy measures and Log reporting
For this version, we've checked two accuracy measurements. In each training/testing iteration (for a mini-batch), we report "true prediction ratio" for each Feature columns and take median accross three Feature categories (e.g. Accetylation, RNA-seq, and TFs).we also report AUC scores in a similar manner. In addition, progress logs are reported in a log file within the current data directory and project name (e.g. deepbrainStaticConvnet) subdirectory.
