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
Data generated from preprocessing steps contains both input DNA seqeunce and corresponding label data for all chromosomes combined in a single file. For training DeepBrain models with this data, we need to seperate them (value-label) and convert them to numpy ndarrays. Moreover, we wanted to leave one chromosome data out while training and test the model with that data. Hence, we have to do train-validation split on both the value and label data based on a single chromosome. This process runs on raijin (for details, see the pbs script: [DL_input_TrainValid_Split_and_Numpy.sh](https://github.com/Akmazad/deepBrain/blob/master/Training%202/pbs%20scripts/DL_input_TrainValid_Split_and_Numpy.sh))

### Training and Validation
 
