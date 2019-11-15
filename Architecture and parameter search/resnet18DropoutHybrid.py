import argparse
import torch
import torch.nn as nn
import math
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
import torch
import torch.nn as nn
import torch.utils.data as tdata
import torch.nn.functional as tfunc
import torch.optim as topti
import logging
import os
import torch.cuda
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix

from functools import partial
from dataclasses import dataclass
from collections import OrderedDict
from torchsummary import summary

import torch.tensor

class config():
  def __init__(self):
    self.name = "hybrid1"
    self.DataDir = "/srv/scratch/z5114185/vafaee/deepBrain/data/"
#     self.TrainingDataFile = "tempTrain5_trainingData_TF_filtered_chr5_value.npy"
#     self.TrainingLabelFile = "tempTrain5_trainingData_TF_filtered_chr5_label_all.npy"
#     self.TestingDataFile = "tempTrain5_validationData_TF_filtered_chr5_value.npy"
#     self.TestingLabelFile = "tempTrain5_validationData_TF_filtered_chr5_label_all.npy"

    self.TrainingDataFile = "HumanFC_OR_EpiMap_AND_ENCODE_tf_specific.bin.Seq_Labels_trainingData_chr1_value.npy"
    self.TrainingLabelFile = "HumanFC_OR_EpiMap_AND_ENCODE_tf_specific.bin.Seq_Labels_trainingData_chr1_label_all.npy"
    self.TestingDataFile = "HumanFC_OR_EpiMap_AND_ENCODE_tf_specific.bin.Seq_Labels_validationData_chr1_value.npy"
    self.TestingLabelFile = "HumanFC_OR_EpiMap_AND_ENCODE_tf_specific.bin.Seq_Labels_validationData_chr1_label_all.npy"
    
    self.w_lr = 1e-2
    self.w_lr_min = 8e-7
    self.w_momentum = 0.9
    self.w_weight_decay = 5e-7
    self.l1_sparsity = 1e-8
    # self.l1_sparsity = 0
    self.print_freq = 5
    self.BATCH_SIZE = 256
    self.val_batch_size = 256
    self.seed = 0
    self.workers = 4 
    self.world_size = -1
    self.rank = 0
    self.dist_url = 'env://'
    self.dist_backend = 'nccl'
    self.gpu = None
    self.multiprocessing_distributed = True
    self.nEpochs = 12
    self.start_epoch = 0

    # architecture-related parameters
    # 3 conv-layers, 1 fully connected layers (See DeepSEA paper)
    self.CONV1_INPUT_CHANNELS = 4
    self.CONV1_OUTPUT_CHANNELS = 320
    self.CONV2_OUTPUT_CHANNELS = 480
    self.CONV3_OUTPUT_CHANNELS = 960
    self.KERNEL_SIZE = 8
    self.POOLING_TH = 4
    self.DROPOUT_l1 = 0.2
    self.DROPOUT_l2 = 0.2
    self.DROPOUT_l3 = 0.5
    self.NUM_OUTPUTS = 566
    self.SEQ_LEN = 1000

    
best_acc1 = 0

# Mohamed's models

'''
	Used for Resnet[18|34] 
'''
class BasicBlock(nn.Module):
    expansion = 1

    def __init__(self, in_planes, planes, stride=1):
        super(BasicBlock, self).__init__()
        self.conv1 = nn.Conv1d(in_planes, planes, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm1d(planes)
        self.conv2 = nn.Conv1d(planes, planes, kernel_size=3, stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm1d(planes)

        self.shortcut = nn.Sequential()
        if stride != 1 or in_planes != self.expansion*planes:
            self.shortcut = nn.Sequential(
                nn.Conv1d(in_planes, self.expansion*planes, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm1d(self.expansion*planes)
            )

    def forward(self, x):
        out = F.relu(self.bn1(self.conv1(x)))
        out = self.bn1(self.conv2(out))
        out += self.shortcut(x)
        out = F.relu(out)
        return out

'''
	Used for Resnet architectures greater than 34 layers
'''
class Bottleneck(nn.Module):
    expansion = 4

    def __init__(self, in_planes, planes, stride=1):
        super(Bottleneck, self).__init__()
        self.conv1 = nn.Conv1d(in_planes, planes, kernel_size=1, bias=False)
        self.bn1 = nn.BatchNorm1d(planes)
        self.conv2 = nn.Conv1d(planes, planes, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn2 = nn.BatchNorm1d(planes)
        self.conv3 = nn.Conv1d(planes, self.expansion*planes, kernel_size=1, bias=False)
        self.bn3 = nn.BatchNorm1d(self.expansion*planes)

        self.shortcut = nn.Sequential()
        if stride != 1 or in_planes != self.expansion*planes:
            self.shortcut = nn.Sequential(
                nn.Conv1d(in_planes, self.expansion*planes, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm1d(self.expansion*planes)
            )

    def forward(self, x):
        out = F.relu(self.bn1(self.conv1(x)))
        out = F.relu(self.bn2(self.conv2(out)))
        out = self.bn3(self.conv3(out))
        out += self.shortcut(x)
        out = F.relu(out)
        return out

'''
	Resnet architecture foundation
'''
class ResNet(nn.Module):
    def __init__(self, block, num_blocks, num_classes=566):
        super(ResNet, self).__init__()
        self.in_planes = 64

        self.conv1 = nn.Conv1d(4, 64, kernel_size=8, stride=1, padding=1, bias=False)
        self.bn1 = nn.BatchNorm1d(64)

        self.layer1 = self._make_layer(block, 64, num_blocks[0], stride=1)
        self.layer2 = self._make_layer(block, 128, num_blocks[1], stride=2)
        self.layer3 = self._make_layer(block, 256, num_blocks[2], stride=2)
        self.layer4 = self._make_layer(block, 512, num_blocks[3], stride=2)
        self.linear = nn.Sequential(nn.Linear(512*block.expansion*31, num_classes)
                                    ,nn.Sigmoid())

    def _make_layer(self, block, planes, num_blocks, stride):
        strides = [stride] + [1]*(num_blocks-1)
        layers = []
        for stride in strides:
            layers.append(block(self.in_planes, planes, stride))
            self.in_planes = planes * block.expansion
        return nn.Sequential(*layers)

    def forward(self, x):
        out = F.relu(self.bn1(self.conv1(x)))
        out = self.layer1(out)
        out = F.dropout(out, p=0.5)
        out = self.layer2(out)
        out = F.dropout(out, p=0.5)
        out = self.layer3(out)
        out = F.dropout(out, p=0.5)
        out = self.layer4(out)
        out = F.dropout(out, p=0.5)
        out = F.avg_pool1d(out, 4)
        out = out.view(out.size(0), -1)
        out = self.linear(out)
        return out


def ResNet18():
    return ResNet(BasicBlock, [2,2,1,1])

def ResNet34():
    return ResNet(BasicBlock, [3,4,6,3])

def ResNet50():
    return ResNet(Bottleneck, [3,4,6,3])

def ResNet101():
    return ResNet(Bottleneck, [3,4,23,3])

def ResNet152():
    return ResNet(Bottleneck, [3,8,36,3])



class ResNet18_2(nn.Module):
    """
    A deeper DeepSEA model architecture.
    Parameters
    ----------
    sequence_length : int
        The length of the sequences on which the model trains and and makes
        predictions.
    n_targets : int
        The number of targets (classes) to predict.
    Attributes
    ----------
    conv_net : torch.nn.Sequential
        The convolutional neural network component of the model.
    classifier : torch.nn.Sequential
        The linear classifier and sigmoid transformation components of the
        model.
    """

    def __init__(self, sequence_length, n_targets):
        super(ResNet18_2, self).__init__()
        conv_kernel_size = 8
        pool_kernel_size = 4

        self.conv_net = nn.Sequential(
            nn.Conv1d(4, 320, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(320, 320, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.BatchNorm1d(320),

            nn.Conv1d(320, 480, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(480, 480, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.BatchNorm1d(480),
            nn.Dropout(p=0.3),

            nn.Conv1d(480, 960, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(960, 960, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(960),
            nn.Dropout(p=0.3))

        reduce_by = 2 * (conv_kernel_size - 1)
        pool_kernel_size = float(pool_kernel_size)
        self._n_channels = int(
            np.floor(
                (np.floor(
                    (sequence_length - reduce_by) / pool_kernel_size)
                 - reduce_by) / pool_kernel_size)
            - reduce_by)
        self.classifier = nn.Sequential(
            nn.Linear(960 * self._n_channels, n_targets),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(n_targets),
            nn.Linear(n_targets, n_targets),
            nn.Sigmoid())

    def forward(self, x):
        """
        Forward propagation of a batch.
        """
        out = self.conv_net(x)
        reshape_out = out.view(out.size(0), 960 * self._n_channels)
        predict = self.classifier(reshape_out)
        return predict


class DeeperDeepSEA(nn.Module):
    """
    A deeper DeepSEA model architecture.
    Parameters
    ----------
    sequence_length : int
        The length of the sequences on which the model trains and and makes
        predictions.
    n_targets : int
        The number of targets (classes) to predict.
    Attributes
    ----------
    conv_net : torch.nn.Sequential
        The convolutional neural network component of the model.
    classifier : torch.nn.Sequential
        The linear classifier and sigmoid transformation components of the
        model.
    """

    def __init__(self, sequence_length, n_targets):
        super(DeeperDeepSEA, self).__init__()
        conv_kernel_size = 8
        pool_kernel_size = 4

        self.conv_net = nn.Sequential(
            nn.Conv1d(4, 320, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(320, 320, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.BatchNorm1d(320),

            nn.Conv1d(320, 480, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(480, 480, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.BatchNorm1d(480),
            nn.Dropout(p=0.3),

            nn.Conv1d(480, 960, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(960, 960, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(960),
            nn.Dropout(p=0.3))

        reduce_by = 2 * (conv_kernel_size - 1)
        pool_kernel_size = float(pool_kernel_size)
        self._n_channels = int(
            np.floor(
                (np.floor(
                    (sequence_length - reduce_by) / pool_kernel_size)
                 - reduce_by) / pool_kernel_size)
            - reduce_by)
        self.classifier = nn.Sequential(
            nn.Linear(960 * self._n_channels, n_targets),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(n_targets),
            nn.Linear(n_targets, n_targets),
            nn.Sigmoid())

    def forward(self, x):
        out = self.conv_net(x)
        reshape_out = out.view(out.size(0), 960 * self._n_channels)
        predict = self.classifier(reshape_out)
        return predict


class LoadDataset(tdata.Dataset):
    def __init__(self, args, dataPath, dataFile, labelFile):
        # Load data from files.
        # self.inputs = np.memmap(dataPath + dataFile, mode="r").reshape(-1, args.CONV1_INPUT_CHANNELS, args.SEQ_LEN)
        # self.labels = np.memmap(dataPath + labelFile, mode="r").reshape(-1, args.NUM_OUTPUTS)

        self.inputs = np.load(dataPath + dataFile)
        self.labels = np.load(dataPath + labelFile)

        self.length = len(self.labels)

    def __getitem__(self, index):
        # Return a single input/label pair from the dataset.
        inputSample = np.array(self.inputs[index], dtype=np.float32)
        labelSample = np.array(self.labels[index], dtype=np.float32)
        sample = (inputSample, labelSample)
        return sample

    def __len__(self):
        return self.length




def getCustomAccuracy2(predicted, target, args):
    n_digits = 3  # to have something like 0.499 = 0.5
    _predicted = torch.round(predicted * 10 ** n_digits) / (10 ** n_digits)
    __predicted = torch.round(_predicted)

    N = predicted.size(0)
    custom_accuracy = np.zeros(args.NUM_OUTPUTS, dtype=np.float)
    for i in range(args.NUM_OUTPUTS):
        truePred = torch.sum(torch.eq(__predicted[:, i], target[:, i])).item()
        custom_accuracy[i] = truePred / N

    return np.median(custom_accuracy[:288]), np.median(custom_accuracy[288:438]), np.median(custom_accuracy[438:])


def getAUCscore(predicted, target, args, logger):
    # n_digits = 3
    # _predicted = torch.round(predicted * 10 ** n_digits) / (10 ** n_digits)
    # __predicted = torch.round(_predicted)

    __predicted = predicted.detach()
    # __predicted = __predicted.detach()
    _target = target.detach()

    aucs = np.empty(args.NUM_OUTPUTS, dtype=float)
    aucs[:] = np.nan

    for i in range(args.NUM_OUTPUTS):
        try:
            auc = roc_auc_score(_target.cpu().numpy()[:, i], __predicted.cpu().numpy()[:, i])
            aucs[i] = auc
        except ValueError as e:
            pass
            # logger.info("NA (No positive (i.e. signal) in Test region)")

    return np.nanmedian(aucs[:288]), np.nanmedian(aucs[288:438]), np.nanmedian(aucs[438:])

def getTPRandTNR(predicted, target, args, logger):
    n_digits = 3
    _predicted = torch.round(predicted * 10 ** n_digits) / (10 ** n_digits)
    __predicted = torch.round(_predicted)

    P = __predicted.cpu()
    T = target.cpu()


    N = T.shape
    #t_ones = torch.sum(torch.eq(torch.ones(N, dtype=torch.float32),T)).item()
    #t_zeros = torch.sum(torch.eq(torch.zeros(N, dtype=torch.float32),T)).item()
    #tpr = torch.sum(torch.eq(T,P) * torch.eq(torch.ones(N, dtype=torch.float32),T)).item()/t_ones
    #tnr = torch.sum(torch.eq(T,P) * torch.eq(torch.zeros(N, dtype=torch.float32),T)).item()/t_zeros
    #return tpr, tnr
    t_ones_H_bool = torch.eq(torch.ones(T.shape[0], 288, dtype=torch.float32),T[:,0:288])
    t_zeros_H_bool = torch.eq(torch.zeros(T.shape[0], 288, dtype=torch.float32),T[:,0:288])
    t_ones_H = torch.sum(t_ones_H_bool).item()
    t_zeros_H = torch.sum(t_zeros_H_bool).item()
    tpr_H = torch.sum(torch.eq(T[:,0:288],P[:,0:288]) * t_ones_H_bool).item()/t_ones_H
    tnr_H = torch.sum(torch.eq(T[:,0:288],P[:,0:288]) * t_zeros_H_bool).item()/t_zeros_H
    
    t_ones_E_bool = torch.eq(torch.ones(T.shape[0], 150, dtype=torch.float32),T[:,288:438])
    t_zeros_E_bool = torch.eq(torch.zeros(T.shape[0], 150, dtype=torch.float32),T[:,288:438])
    t_ones_E = torch.sum(t_ones_E_bool).item()
    t_zeros_E = torch.sum(t_zeros_E_bool).item()
    tpr_E = torch.sum(torch.eq(T[:,288:438],P[:,288:438]) * t_ones_E_bool).item()/t_ones_E
    tnr_E = torch.sum(torch.eq(T[:,288:438],P[:,288:438]) * t_zeros_E_bool).item()/t_zeros_E
  
    t_ones_T_bool = torch.eq(torch.ones(T.shape[0], 128, dtype=torch.float32),T[:,438:])
    t_zeros_T_bool = torch.eq(torch.zeros(T.shape[0], 128, dtype=torch.float32),T[:,438:])
    t_ones_T = torch.sum(t_ones_T_bool).item()
    t_zeros_T = torch.sum(t_zeros_T_bool).item()
    tpr_T = torch.sum(torch.eq(T[:,438:],P[:,438:]) * t_ones_T_bool).item()/t_ones_T
    tnr_T = torch.sum(torch.eq(T[:,438:],P[:,438:]) * t_zeros_T_bool).item()/t_zeros_T
    
    t_ones_bool = torch.cat((t_ones_H_bool,t_ones_E_bool,t_ones_T_bool),1)
    t_zeros_bool = torch.cat((t_zeros_H_bool,t_zeros_E_bool,t_zeros_T_bool),1)
    t_ones = torch.sum(t_ones_bool).item()
    t_zeros = torch.sum(t_zeros_bool).item()
    tpr_all = torch.sum(torch.eq(T,P) * t_ones_bool).item()/t_ones
    tnr_all = torch.sum(torch.eq(T,P) * t_zeros_bool).item()/t_zeros


    return [tpr_all, tnr_all, tpr_H, tnr_H, tpr_E, tnr_E, tpr_T, tnr_T]


def getAccuracyScore(predicted, target, args, logger):
    n_digits = 3
    _predicted = torch.round(predicted * 10**n_digits) / (10**n_digits)
    __predicted = torch.round(_predicted)

    # _predicted = torch.round(predicted).detach()
    __predicted = __predicted.detach()
    _target = target.detach()

    accs = np.zeros(args.NUM_OUTPUTS, dtype=np.float)
    for i in range(args.NUM_OUTPUTS):
        try:
            acc = accuracy_score(_target.cpu().numpy()[:, i], __predicted.cpu().numpy()[:, i])
            accs[i] = acc
        except ValueError as e:
            pass
            # logger.info("NA (No positive (i.e. signal) in Test region)")

    return np.median(accs[:288]), np.median(accs[288:438]), np.median(accs[438:])


def getF1Score(predicted, target, args, logger):
    n_digits = 3
    _predicted = torch.round(predicted * 10 ** n_digits) / (10 ** n_digits)
    __predicted = torch.round(_predicted)

    #     _predicted = torch.round(predicted).detach()
    __predicted = __predicted.detach()
    _target = target.detach()

    N = target.size(0)

    f1s = np.zeros(args.NUM_OUTPUTS, dtype=np.float)
    for i in range(args.NUM_OUTPUTS):

        nOnes = torch.sum(torch.eq(torch.ones(N), target.cpu()[:, i])).item()
        nZeros = torch.sum(torch.eq(torch.zeros(N), target.cpu()[:, i])).item()
        if (nOnes != N & nZeros != N):
            f1 = f1_score(_target.cpu().numpy()[:, i], __predicted.cpu().numpy()[:, i])
            f1s[i] = f1
        else:
            pass
            # logger.info("NA (No positive (i.e. signal) in Test region)")

    return np.median(f1s[:288]), np.median(f1s[288:438]), np.median(f1s[438:])


def get_logger(file_path):
    """ Make python logger """
    # [!] Since tensorboardX use default logger (e.g. logging.info()), we should use custom logger
    logger = logging.getLogger('db2')
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


def find_perc_uncertainty(output, low, high):
    n_digits = 3
    output = torch.round(output * 10**n_digits) / (10**n_digits)
    output = output.detach().cpu().numpy()
    return np.sum(np.logical_and(output>=low, output<=high))/output.size    # returns the proportion of tensor elements are within the range

def train(train_loader, model, criterion, optimizer, epoch, args, logger, device):

    # switch to train mode
    model.train()
    # perc_uncertainty = 0.0

    for i, (input, target) in enumerate(train_loader):
      
        input = input.cuda(device, non_blocking=True)
        target = target.cuda(device, non_blocking=True)
        
        # compute output
        # output = model(input, args)
        output = model(input)
        loss = criterion(output, target.squeeze(1))
        temp_target = target
        temp_pred = output
        nOnes_target = torch.sum(temp_target.detach()).item()
        temp_pred = torch.round(temp_pred * 10 ** 3) / (10 ** 3)
        temp_pred = torch.round(temp_pred)
        nOnes_predicted = torch.sum(temp_pred.detach()).item()
        tprs_tnrs = getTPRandTNR(output, target, args, logger)


        # measure accuracy and record loss
        # acc = getCustomAccuracy(output, target, args)
        perc_uncertainty = find_perc_uncertainty(output, 0.4, 0.6)
        custom_accuracy = getCustomAccuracy2(output, target, args)
        aucs = getAUCscore(output, target, args, logger)
        accs = getAccuracyScore(output, target, args, logger)
        # f1s = getF1Score(output, target, args, logger)

        # compute gradient
        optimizer.zero_grad()

        # add l1 Sparsity
        l1 = 0
        for p in model.parameters():
            l1 = l1 + p.abs().sum()
        loss = loss + args.l1_sparsity * l1

        # do SGD step
        loss.backward()
        optimizer.step()

        if i % args.print_freq == 0 or i == len(train_loader) - 1:
            #           logger.info("TRAINING: Epoch: %d, Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[HumanFC:%.3f, EpiMap:%.3f, TFs:%.3f], roc[HumanFC:%.3f, EpiMap:%.3f, TFs:%.3f]"
            #                       % (epoch + 1, i, len(train_loader) - 1,
            #                          loss, perc_uncertainty, custom_accuracy[0],
            #                          custom_accuracy[1], custom_accuracy[2],
            #                          aucs[0], aucs[1], aucs[2]))
            logger.info(
                # "TRAINING: Epoch: %d, Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f], roc[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], accuracy[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], f1[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ]"
                "TRAINING: Epoch: %d, Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, nOnes_target: %d, nOnes_pred: %d, TPR[all: %.3f, HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f], TNR[all: %.3f, HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f], roc[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], accuracy[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ]"
                % (epoch + 1, i, len(train_loader) - 1, 
                   loss, perc_uncertainty, nOnes_target, nOnes_predicted, 
                   tprs_tnrs[0], tprs_tnrs[2],tprs_tnrs[4], tprs_tnrs[6],
                   tprs_tnrs[1], tprs_tnrs[3],tprs_tnrs[5], tprs_tnrs[7],
                   aucs[0], aucs[1], aucs[2],
                   accs[0], accs[1], accs[2]))


def validate(val_loader, model, criterion, args, logger, device):

    # switch to evaluate mode
    model.eval()
    total_ACC, total_RNA, total_TFs = 0, 0, 0

    # losses = []

    with torch.no_grad():
        for i, (input, target) in enumerate(val_loader):
            
            input = input.cuda(device, non_blocking=True)
            target = target.cuda(device, non_blocking=True)
        
            output = model(input)
            loss = criterion(output, target)
            temp_target = target
            temp_pred = output
            nOnes_target = torch.sum(temp_target.detach()).item()
            temp_pred = torch.round(temp_pred * 10 ** 3) / (10 ** 3)
            temp_pred = torch.round(temp_pred)
            nOnes_predicted = torch.sum(temp_pred.detach()).item()
            tprs_tnrs = getTPRandTNR(output, target, args, logger)


            # losses[i] = loss

            # measure accuracy and record loss
            # acc = getCustomAccuracy(output, target, args)
            p = find_perc_uncertainty(output, 0.4, 0.6)
            custom_accuracy = getCustomAccuracy2(output, target, args)
            aucs = getAUCscore(output, target, args, logger)
            accs = getAccuracyScore(output, target, args, logger)
            # f1s = getF1Score(output, target, args, logger)

            total_ACC += np.median(aucs[0])
            total_RNA += np.median(aucs[1])
            total_TFs += np.median(aucs[2])

            if i % args.print_freq == 0 or i == len(val_loader) - 1:
                # progress._print(i)
                # logger.info("batch: %d, loss: %.3f; valid accuracy: custom_accuracy_metric: %.3f, ACC marks: %.3f, RNA-seq: %.3f, TFs: %.3f" % (i+1, loss, acc, aucs[0], aucs[1], aucs[2]))
                #                 logger.info("VALIDATION: Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[HumanFC:%.3f, EpiMap:%.3f, TFs:%.3f], roc[HumanFC: %.3f, EpiMap: %.3f, TFs:%.3f]" % (i, len(val_loader)-1, loss, perc_uncertainty, custom_accuracy[0], custom_accuracy[1], custom_accuracy[2], aucs[0], aucs[1], aucs[2]))
                logger.info(
                    "VALIDATION: Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, nOnes_target: %d, nOnes_pred: %d, TPR[all: %.3f, HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f], TNR[all: %.3f, HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f], roc[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], accuracy[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ]"
                    # "VALIDATION: Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f], roc[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], accuracy[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], f1[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ]"
                    % (i, len(val_loader) - 1,  
                       loss, p, nOnes_target, nOnes_predicted, 
                       tprs_tnrs[0], tprs_tnrs[2],tprs_tnrs[4], tprs_tnrs[6],
                       tprs_tnrs[1], tprs_tnrs[3],tprs_tnrs[5], tprs_tnrs[7],
                       aucs[0], aucs[1], aucs[2],
                       accs[0], accs[1], accs[2]))

    total_ACC /= len(val_loader)
    total_RNA /= len(val_loader)
    total_TFs /= len(val_loader)

    return np.median([total_ACC, total_RNA, total_TFs])
    # return np.median(losses), np.median([total_ACC, total_RNA, total_TFs])


def adjust_learning_rate(optimizer, epoch, args, lr_scheduler):
    """Sets the learning rate to the initial LR decayed by 10 every 30 epochs"""
    # lr = args.w_lr * (0.1 ** (epoch // 30))   # option 1
    lr_scheduler.step()     # option 2
    lr = lr_scheduler.get_lr()[0]

    for param_group in optimizer.param_groups:
        param_group['lr'] = lr


def main():
    # Training settings
    args = config()
    use_cuda = torch.cuda.is_available()
    torch.manual_seed(args.seed)
    device = torch.device("cuda" if use_cuda else "cpu")
    print(device)

    DataPath = args.DataDir
    if not (os.path.exists(os.path.join(DataPath, args.name))):
        os.mkdir(os.path.join(DataPath, args.name))
    f = open(os.path.join(DataPath, args.name)+"/"+"{}.log".format(args.name), "a")
    f.close()
    logger = get_logger(os.path.join(os.path.join(DataPath, args.name), "{}.log".format(args.name)))

    trainDataset = LoadDataset(args, dataPath=DataPath, dataFile=args.TrainingDataFile,
                               labelFile=args.TrainingLabelFile)
    # define sampler
    train_sampler = torch.utils.data.sampler.RandomSampler(trainDataset)

    trainLoader = tdata.DataLoader(trainDataset, batch_size=args.BATCH_SIZE, shuffle=(train_sampler is None),
                                   num_workers=args.workers, pin_memory=True, sampler=train_sampler)

    # Load the testing dataset, and create a data loader to generate a batch. CHANGE THIS ONCE TESTING DATASET IS READY
    valDataset = LoadDataset(args, dataPath=DataPath, dataFile=args.TestingDataFile, labelFile=args.TestingLabelFile)
    valLoader = tdata.DataLoader(dataset=valDataset, batch_size=args.val_batch_size, shuffle=False,
                                 num_workers=args.workers, pin_memory=True)

    # build the model and criterion
    # model = BuildModel(args).to(device)
    # model = DeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS).to(device)
    # model = DeeperDeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS).to(device)
    # model = get_model(load_weights=False)
    # model = ResNet18_2(args.SEQ_LEN, args.NUM_OUTPUTS).to(device)
    model = ResNet18().to(device)
    # Add a sigmoid activation function to the output.  Use a binary cross entropy
    criterion = nn.BCELoss().to(device)
    # criterion = nn.CrossEntropyLoss().to(device)

    # optimiser = topti.Adam(model.parameters(), lr=args.w_lr)  # Minimise the loss using the Adam algorithm.
#     optimiser = torch.optim.Adam(model.parameters(), args.w_lr, betas=(0.5, 0.999),
#                                    weight_decay=args.w_weight_decay)
    # print(model.parameters())
    optimiser = torch.optim.SGD(model.parameters(), args.w_lr, momentum=args.w_momentum, weight_decay=args.w_weight_decay)

    # model = DeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS).to(device)
    # criterion = nn.BCELoss()
    # optimiser = (torch.optim.SGD, {"lr": args.w_lr, "weight_decay": 1e-6, "momentum": 0.9})



    lr_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimiser, args.nEpochs, eta_min=args.w_lr_min)
    # lr_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimiser, args.nEpochs, eta_min=0)
    # scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimiser, 'min')

    best_acc1 = 0.0
    for epoch in range(args.start_epoch, args.nEpochs):

        adjust_learning_rate(optimiser, epoch, args, lr_scheduler)

        # train for one epoch
        train(trainLoader, model, criterion, optimiser, epoch, args, logger, device)

        # evaluate on validation set
        acc1 = validate(valLoader, model, criterion, args, logger, device)

        # scheduler.step(acc1[1])

        # remember best acc@1 and save checkpoint
        is_best = acc1 > best_acc1
        best_acc1 = max(acc1, best_acc1)

    # if (args.save_model):
    #     torch.save(model.state_dict(), "mnist_cnn.pt")
    torch.save(model.state_dict(), "deepBrain2_cnn.pt")


if __name__ == '__main__':
    main()
