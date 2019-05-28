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
| [DL_model_test.py](https://github.com/Akmazad/deepBrain/blob/master/Training%202/DL_model_test.py) | DeepBrain with static Convnet (e.g. [DeepSEA](https://github.com/FunctionLab/selene/blob/master/models/deepsea.py) and [DeeperDeepSEA](https://github.com/FunctionLab/selene/blob/master/tutorials/quickstart_training/deeperdeepsea.py)) models |
| [DL_model_test_DARTS.py](https://github.com/Akmazad/deepBrain/blob/master/Training%202/DL_model_test_DARTS.py) | [DARTs](https://github.com/quark0/darts) (Differentiable Architecture Search) implementation of the "DL_model_test.py" * | 
