""" Utilities """
import os
import logging
import shutil
import torch
import torchvision.datasets as dset
import numpy as np
import preproc

import torch.utils.data as tdata

# Class for reading training/testing dataset files.
class deepBrainData(tdata.Dataset):
    def __init__(self, dataFile, labelFile, dataPath):
        # Load data from files.
        self.inputs = np.memmap(dataPath + dataFile, mode = "r").reshape(-1, 4, 1000)
        self.labels = np.memmap(dataPath + labelFile, mode = "r").reshape(-1, 131)

        self.length = len(self.labels)

    def __getitem__(self, index):
        # Return a single input/label pair from the dataset.
        inputSample = np.array(self.inputs[index], dtype = np.float32)
        labelSample = np.array(self.labels[index], dtype = np.float32)
        sample = (inputSample, labelSample)

        return sample

    def __len__(self):

        return self.length

def get_data(train_data, train_label, data_path, logger, cutout_length, validation):
    #logger.info('train_label size: {}'.format(train_label.size(1)))
    # n_classes = train_label.size(1)
    n_classes = 131
    trn_data = deepBrainData(train_data, train_label, data_path)
    shape = trn_data.inputs.shape
    input_channels = shape[1]
    input_size = shape[2]
    # logger.info('input size: {}'.format(input_size))
    # logger.info('n_classes: {}'.format(train_label.size(1)))
    # logger.info('input channels: {}'.format(input_size))



    ret = [input_size, input_channels, n_classes, trn_data]

    return ret


def get_logger(file_path):
    """ Make python logger """
    # [!] Since tensorboardX use default logger (e.g. logging.info()), we should use custom logger
    logger = logging.getLogger('darts')
    log_format = '%(asctime)s | %(message)s'
    formatter = logging.Formatter(log_format, datefmt='%m/%d %I:%M:%S %p')
    file_handler = logging.FileHandler(file_path)
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(logging.INFO)

    return logger


def param_size(model):
    """ Compute parameter size in MB """
    n_params = sum(
        np.prod(v.size()) for k, v in model.named_parameters() if not k.startswith('aux_head'))
    return n_params / 1024. / 1024.


class AverageMeter():
    """ Computes and stores the average and current value """
    def __init__(self):
        self.reset()

    def reset(self):
        """ Reset all statistics """
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        """ Update statistics """
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count


class SumMeter():
    """ Computes and stores the average and current value """
    def __init__(self):
        self.reset()

    def reset(self):
        """ Reset all statistics """
        self.truePos = 0
        self.trueNeg = 0
        self.falsePos = 0
        self.falseNeg = 0

    def update(self, truePos, trueNeg, falsePos, falseNeg):
        """ Update statistics """
        self.truePos += truePos
        self.trueNeg += trueNeg
        self.falsePos += falsePos
        self.falseNeg += falseNeg

    def accuracy(self):
        total = self.truePos + self.trueNeg + self.falsePos + self.falseNeg
        if total==0:
            return 0
        return (self.truePos + self.trueNeg) / total

    def MCC(self):
        numerator = self.truePos * self.trueNeg - self.falsePos * self.falseNeg
        denominator = ((self.truePos + self.falsePos) * (self.truePos + self.falseNeg) * (self.trueNeg + self.falsePos) * (self.trueNeg + self.falseNeg)) ** 0.5

        with np.errstate(divide = "ignore", invalid = "ignore"):
            return np.divide(numerator, denominator)


def accuracy(output, target):
    # """ Computes the precision@k for the specified values of k """
    # maxk = max(topk)
    # batch_size = target.size(0)

    # _, pred = output.topk(maxk, 1, True, True)
    # pred = pred.t()
    # # one-hot case
    # if target.ndimension() > 1:
    #     target = target.max(1)[1]

    # correct = pred.eq(target.view(1, -1).expand_as(pred))

    # res = []
    # for k in topk:
    #     correct_k = correct[:k].view(-1).float().sum(0)
    #     res.append(correct_k.mul_(1.0 / batch_size))

    # return res

    predicted = torch.round(torch.sigmoid(output))

    truePos = torch.sum(target * predicted).item()
    trueNeg = torch.sum((1 - target) * (1 - predicted)).item()
    falsePos = torch.sum((1 - target) * predicted).item()
    falseNeg = torch.sum(target * (1 - predicted)).item()

    return [truePos, trueNeg, falsePos, falseNeg]


def save_checkpoint(state, ckpt_dir, is_best=False):
    filename = os.path.join(ckpt_dir, 'checkpoint.pth.tar')
    torch.save(state, filename)
    if is_best:
        best_filename = os.path.join(ckpt_dir, 'best.pth.tar')
        shutil.copyfile(filename, best_filename)
