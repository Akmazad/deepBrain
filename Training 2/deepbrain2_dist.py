import argparse
import os
import random
import shutil
import time
import warnings
import sys
import math
import numpy as np
import logging
from sklearn.metrics import roc_auc_score

import torch
import torch.nn as nn
import torch.nn.parallel
import torch.optim as topti
import torch.nn.functional as tfunc

import torch.backends.cudnn as cudnn
import torch.distributed as dist
import torch.multiprocessing as mp
import torch.utils.data as tdata
import torch.utils.data.distributed

from torch.cuda.comm import broadcast_coalesced
from torch.cuda import nccl

if dist.is_available():
    from torch.distributed.distributed_c10d import _get_default_group


from torch.cuda._utils import _get_device_index

# configuration parameters
parser = argparse.ArgumentParser("deepBrain2")
parser.add_argument('--name', required=True)
parser.add_argument('--DataDir', required=True)
parser.add_argument('--TrainingDataFile', required=True)
parser.add_argument('--TrainingLabelFile', required=True)
parser.add_argument('--TestingDataFile', required=True)
parser.add_argument('--TestingLabelFile', required=True)
# parser.add_argument('--w_lr', type=float, default=0.01, help='lr for weights')
parser.add_argument('--w_lr', type=float, default=1e-2, help='lr for weights')
parser.add_argument('--w_lr_min', type=float, default=8e-7, help='minimum lr for weights')
parser.add_argument('--w_momentum', type=float, default=0.9)
# parser.add_argument('--w_weight_decay', type=float, default=3e-4)
parser.add_argument('--w_weight_decay', type=float, default=5e-7)
parser.add_argument('--w_grad_clip', type=float, default=5.,
                    help='gradient clipping for weights')
parser.add_argument('--print_freq', type=int, default=50, help='print frequency')
parser.add_argument('--init_channels', type=int, default=16)
parser.add_argument('--layers', type=int, default=8)
parser.add_argument('--BATCH_SIZE', type=int, default=256, metavar='N',
                    help='mini-batch size, this is the total batch size of all GPUs on the current '
                         'node when using Data Parallel or Distributed Data Parallel')
parser.add_argument('--seed', type=int, default=None, help='seed for initializing training')
parser.add_argument('--workers', type=int, default=4, help='# of workers')
parser.add_argument('--alpha_lr', type=float, default=3e-4, help='lr for alpha')
parser.add_argument('--alpha_weight_decay', type=float, default=1e-4, help='weight decay for alpha')
parser.add_argument('--world-size', default=-1, type=int, help='number of nodes for distributed training')
parser.add_argument('--rank', default=0, type=int, help='node rank for distributed training')
parser.add_argument('--dist-url', default='env://', type=str, help='url used to set up distributed training')
parser.add_argument('--dist-backend', default='gloo', type=str, help='distributed backend')
parser.add_argument('--gpu', default=None, type=int, help='GPU id to use.')
parser.add_argument('--multiprocessing-distributed', default=True, action='store_true',
                    help='Use multi-processing distributed training to launch N processes per node, '
                         'which has N GPUs. This is the fastest way to use PyTorch for either single '
                         'node or multi node data parallel training')
parser.add_argument('--pretrained', dest='pretrained', action='store_true', help='use pre-trained model')
parser.add_argument('--nEpochs', default=5, type=int, metavar='N', help='number of total epochs to run')
parser.add_argument('--start-epoch', default=0, type=int, metavar='N', help='manual epoch number (useful on restarts)')
parser.add_argument('--resume', default='', type=str, metavar='PATH', help='path to latest checkpoint (default: none)')

# architecture-related parameters
# 3 conv-layers, 1 fully connected layers (See DeepSEA paper)
parser.add_argument('--CONV1_INPUT_CHANNELS', type=int, default=4)
parser.add_argument('--CONV1_OUTPUT_CHANNELS', type=int, default=320)
parser.add_argument('--CONV2_OUTPUT_CHANNELS', type=int, default=480)
parser.add_argument('--CONV3_OUTPUT_CHANNELS', type=int, default=960)
parser.add_argument('--KERNEL_SIZE', type=int, default=8)
parser.add_argument('--POOLING_TH', type=int, default=4)
parser.add_argument('--DROPOUT_l1', type=float, default=0.2)
parser.add_argument('--DROPOUT_l2', type=float, default=0.2)
parser.add_argument('--DROPOUT_l3', type=float, default=0.5)
parser.add_argument('--NUM_OUTPUTS', type=int, default=131)
parser.add_argument('--SEQ_LEN', type=int, default=1000)

best_acc1 = 0


# MODEL zone

# Own model (DeepSEA like)
class BuildModel(nn.Module):
    def __init__(self, args):
        super(BuildModel, self).__init__()

        # Create and initialise weights and biases for the layers.
        # regularization parameters could be introduced later (see deepsea/deepsea_model.py)
        self.conv_layer1 = nn.Conv1d(args.CONV1_INPUT_CHANNELS, args.CONV1_OUTPUT_CHANNELS, args.KERNEL_SIZE)
        self.conv_layer2 = nn.Conv1d(args.CONV1_OUTPUT_CHANNELS, args.CONV2_OUTPUT_CHANNELS, args.KERNEL_SIZE)
        self.conv_layer3 = nn.Conv1d(args.CONV2_OUTPUT_CHANNELS, args.CONV3_OUTPUT_CHANNELS, args.KERNEL_SIZE)

        nChannel = math.floor((args.SEQ_LEN - (args.KERNEL_SIZE - 1)) / args.POOLING_TH)
        nChannel = math.floor((nChannel - (args.KERNEL_SIZE - 1)) / args.POOLING_TH)
        nChannel = math.floor((nChannel - (args.KERNEL_SIZE - 1)) / args.POOLING_TH)
        self.fc1 = nn.Linear(args.CONV3_OUTPUT_CHANNELS * nChannel, args.NUM_OUTPUTS)

    def forward(self, x, args):
        # Create the forward pass through the network.
        x = self.conv_layer1(x)
        x = tfunc.leaky_relu(x, 0.01)
        x = tfunc.max_pool1d(x, args.POOLING_TH)  # downsample by half (i.e. if parameter=4, then by quarter)
        x = tfunc.dropout(x, args.DROPOUT_l1)

        x = self.conv_layer2(x)
        x = tfunc.leaky_relu(x, 0.01)
        x = tfunc.max_pool1d(x, args.POOLING_TH)
        x = tfunc.dropout(x, args.DROPOUT_l2)

        x = self.conv_layer3(x)
        x = tfunc.leaky_relu(x, 0.1)
        x = tfunc.max_pool1d(x, args.POOLING_TH)
        x = tfunc.dropout(x, args.DROPOUT_l3)

        # for the fully connected layer
        x = x.view(x.shape[0], -1)  # Flatten tensor.
        x = self.fc1(x)
        x = tfunc.leaky_relu(x, 0.01)
        #         x = tfunc.dropout(x, 0.2)
        x = torch.sigmoid(x)

        return x

# DeepSEA (copied from "FunctionLab/selene")
class DeepSEA(nn.Module):
    def __init__(self, sequence_length, n_genomic_features):
        """
        Parameters
        ----------
        sequence_length : int
        n_genomic_features : int
        """
        super(DeepSEA, self).__init__()
        conv_kernel_size = 8
        pool_kernel_size = 4

        self.conv_net = nn.Sequential(
            nn.Conv1d(4, 320, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.Dropout(p=0.2),

            nn.Conv1d(320, 480, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.Dropout(p=0.2),

            nn.Conv1d(480, 960, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.5))

        reduce_by = conv_kernel_size - 1
        pool_kernel_size = float(pool_kernel_size)
        self.n_channels = int(
            np.floor(
                (np.floor(
                    (sequence_length - reduce_by) / pool_kernel_size)
                 - reduce_by) / pool_kernel_size)
            - reduce_by)
        self.classifier = nn.Sequential(
            nn.Linear(960 * self.n_channels, n_genomic_features),
            nn.ReLU(inplace=True),
            nn.Linear(n_genomic_features, n_genomic_features),
            nn.Sigmoid())

    def forward(self, x):
        """Forward propagation of a batch.
        """
        out = self.conv_net(x)
        reshape_out = out.view(out.size(0), 960 * self.n_channels)
        predict = self.classifier(reshape_out)
        return predict

# deeperDeepSEA (copied from "FunctionLab/selene")
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
            nn.Dropout(p=0.2),

            nn.Conv1d(480, 960, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(960, 960, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(960),
            nn.Dropout(p=0.2))

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
# Model Zone -- Ends

# Loads data into numpy ndarray
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


# function
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


"""
  This function returns the true prediction ratio for each features
  
  Parameters
  ----------
  'predicted': 1D tensor (for a batch)
  'target': 1D tensor (for a batch)
  'args': argument vector (for a batch)
  
  Return
  ------
  An array of median accuracies accross each feature categories (Accetylation, RNS-seq, and TFs
"""

def getCustomAccuracy2(predicted, target, args):
    # predicted = torch.round(torch.sigmoid(predicted))
    # predicted = torch.round(predicted)
    n_digits = 3    # to have something like 0.499 = 0.5
    _predicted = torch.round(predicted * 10 ** n_digits) / (10 ** n_digits)
    __predicted = torch.round(_predicted)

    N = predicted.size(0)
    custom_accuracy = np.zeros(args.NUM_OUTPUTS, dtype=np.float)
    for i in range(args.NUM_OUTPUTS):
        truePred = torch.sum(torch.eq(__predicted[:, i], target[:, i])).item()
        custom_accuracy[i] = truePred/N

    return np.median(custom_accuracy[:2]), np.median(custom_accuracy[2]), np.median(custom_accuracy[3:])


"""
  This function returns the AUC scores for each features
  
  Parameters
  ----------
  'predicted': 1D tensor (for a batch)
  'target': 1D tensor (for a batch)
  'args': argument vector (for a batch)
  'logger': an object for keeping progress logs
  
  Return
  ------
  An array of median AUCs accross each feature categories (Accetylation, RNS-seq, and TFs
"""

def getAUCscore(predicted, target, args, logger):
    # n_digits = 3
    # _predicted = torch.round(predicted * 10**n_digits) / (10**n_digits)
    # __predicted = torch.round(_predicted)

    # _predicted = torch.round(predicted).detach()
    __predicted = predicted.detach()
    _target = target.detach()

    aucs = np.zeros(args.NUM_OUTPUTS, dtype=np.float)
    for i in range(args.NUM_OUTPUTS):
        try:
            auc = roc_auc_score(_target.cpu().numpy()[:, i], __predicted.cpu().numpy()[:, i], average='weighted')
            aucs[i] = auc
        except ValueError as e:
            pass
            # logger.info("NA (No positive (i.e. signal) in Test region)")

    return np.median(aucs[:2]), np.median(aucs[2]), np.median(aucs[3:])


# This function is a modified version of this: https://github.com/pytorch/examples/blob/master/imagenet/main.py
def main():
    # distributed code
    try:
        mp.set_start_method('spawn')
    except RuntimeError:
        pass

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        torch.manual_seed(args.seed)
        cudnn.deterministic = True
        warnings.warn('You have chosen to seed training. '
                      'This will turn on the CUDNN deterministic setting, '
                      'which can slow down your training considerably! '
                      'You may see unexpected behavior when restarting '
                      'from checkpoints.')

    if args.gpu is not None:
        warnings.warn('You have chosen a specific GPU. This will completely '
                      'disable data parallelism.')

    if args.dist_url == "env://" and args.world_size == -1:
        args.world_size = int(os.environ["WORLD_SIZE"])

    args.distributed = args.world_size > 1 or args.multiprocessing_distributed

    ngpus_per_node = torch.cuda.device_count()
    try:
        if args.multiprocessing_distributed:
            # Since we have ngpus_per_node processes per node, the total world_size
            # needs to be adjusted accordingly
            args.world_size = ngpus_per_node * args.world_size
            # Use torch.multiprocessing.spawn to launch distributed processes: the
            # main_worker process function
            mp.spawn(main_worker, nprocs=ngpus_per_node, args=(ngpus_per_node, args))
        else:
            # Simply call main_worker function
            main_worker(args.gpu, ngpus_per_node, args)
    except RuntimeError as e:
        print(e)

# Each GPU runs this process worker
def main_worker(gpu, ngpus_per_node, args):
    global best_acc1
    args.gpu = gpu

    DataPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    DataPath = os.path.join(DataPath, args.DataDir)
    logger = get_logger(os.path.join(os.path.join(DataPath, args.name, 'staticConvNet'), "{}.log".format(args.name)))

    if args.gpu is not None:
        print("Use GPU: {} for training".format(args.gpu))

    if args.distributed:
        if args.dist_url == "env://" and args.rank == -1:
            args.rank = int(os.environ["RANK"])
        if args.multiprocessing_distributed:
            # For multiprocessing distributed training, rank needs to be the
            # global rank among all the processes
            args.rank = args.rank * ngpus_per_node + gpu
            
        # intialize the process group in the distributed environment
        dist.init_process_group(backend=args.dist_backend, init_method=args.dist_url,
                                world_size=args.world_size, rank=args.rank)

    # build/Select the model and criterion
    # model = BuildModel(args)
    # model = DeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS)
    model = DeeperDeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS)

    if args.distributed:
        # For multiprocessing distributed, DistributedDataParallel constructor
        # should always set the single device scope, otherwise,
        # DistributedDataParallel will use all available devices.
        if args.gpu is not None:
            torch.cuda.set_device(args.gpu)
            model.cuda(args.gpu)
            # When using a single GPU per process and per
            # DistributedDataParallel, we need to divide the batch size
            # ourselves based on the total number of GPUs we have
            args.BATCH_SIZE = int(args.BATCH_SIZE / ngpus_per_node)
            args.workers = int(args.workers / ngpus_per_node)
            logger.info("args.workers: {}, and args.gpu: {}".format(args.workers, args.gpu))
            model = torch.nn.parallel.DistributedDataParallel(model, device_ids=[args.gpu])
        else:
            model.cuda()
            # DistributedDataParallel will divide and allocate BATCH_SIZE to all
            # available GPUs if device_ids are not set
            model = torch.nn.parallel.DistributedDataParallel(model)
    elif args.gpu is not None:
        torch.cuda.set_device(args.gpu)
        model = model.cuda(args.gpu)
    else:
        # DataParallel will divide and allocate BATCH_SIZE to all available GPUs
        if args.arch.startswith('alexnet') or args.arch.startswith('vgg'):
            model.features = torch.nn.DataParallel(model.features)
            model.cuda()
        else:
            model = torch.nn.DataParallel(model).cuda()

    # Add a sigmoid activation function to the output.  Use a binary cross entropy
    # criterion = nn.BCEWithLogitsLoss()
    criterion = nn.BCELoss().cuda(args.gpu)
    # loss function.
    # optimiser = topti.Adam(model.parameters(), lr=args.w_lr)  # Minimise the loss using the Adam algorithm.
    optimiser = torch.optim.SGD(model.parameters(), args.w_lr, momentum=args.w_momentum,
                              weight_decay=args.w_weight_decay)

    cudnn.benchmark = True

    # Load the training dataset, and create a data loader to generate a batch.
    trainDataset = LoadDataset(args, dataPath=DataPath, dataFile=args.TrainingDataFile, labelFile=args.TrainingLabelFile)
    # define sampler
    if args.distributed:
        train_sampler = torch.utils.data.distributed.DistributedSampler(trainDataset)
    else:
        train_sampler = None

    trainLoader = tdata.DataLoader(trainDataset, batch_size=args.BATCH_SIZE, shuffle=(train_sampler is None),
                                   num_workers=args.workers, pin_memory=True, sampler=train_sampler)

    # Load the testing dataset, and create a data loader to generate a batch. CHANGE THIS ONCE TESTING DATASET IS READY
    valDataset = LoadDataset(args, dataPath=DataPath, dataFile=args.TestingDataFile, labelFile=args.TestingLabelFile)
    valLoader = tdata.DataLoader(dataset=valDataset, batch_size=args.BATCH_SIZE, shuffle=False,
                                   num_workers=args.workers, pin_memory=True)

    lr_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimiser, args.nEpochs, eta_min=args.w_lr_min)
    # lr_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimiser, args.nEpochs, eta_min=0)

    logger.info("Training starts: ")

    for epoch in range(args.start_epoch, args.nEpochs):
        if args.distributed:
            train_sampler.set_epoch(epoch)

        # learning rate adjusting
        adjust_learning_rate(optimiser, epoch, args, lr_scheduler)    # option1

        # train for one epoch
        train(trainLoader, model, criterion, optimiser, epoch, args, logger)

        # evaluate on validation set
        acc1 = validate(valLoader, model, criterion, args, logger)


"""
  This function returns "percentage of uncertain predictions" e.g. range: [0.4, 0.6]
  
  Parameters
  ----------
  'output': a vector (for a batch)
  'low': a vector (for a batch)
  'high': argument vector (for a batch)
  Return
  ------
  The proportion of tensor elements are within the given range
"""
        
def find_perc_uncertainty(output, low, high):
    output = output.detach().cpu().numpy()
    return np.sum(np.logical_and(output >= low, output <= high)) / output.size    # returns the proportion of tensor elements are within the range


def train(train_loader, model, criterion, optimizer, epoch, args, logger):
    batch_time = AverageMeter('Time', ':6.3f')
    data_time = AverageMeter('Data', ':6.3f')
    
    # switch to train mode
    model.train()

    end = time.time()
    for i, (input, target) in enumerate(train_loader):
        # measure data loading time
        data_time.update(time.time() - end)

        if args.gpu is not None:
            input = input.cuda(args.gpu, non_blocking=True)
            target = target.cuda(args.gpu, non_blocking=True)

        # compute output
        # output = model(input, args)
        output = model(input)
        loss = criterion(output, target)

        # measure accuracy and record loss
        perc_uncertainty = find_perc_uncertainty(output, 0.4, 0.6)
        custom_accuracy = getCustomAccuracy2(output, target, args)
        aucs = getAUCscore(output, target, args, logger)

        # compute gradient
        optimizer.zero_grad()

        # add l1 Sparsity
        l1 = 0
        for p in model.parameters():
            l1 = l1 + p.abs().sum()
        loss = loss + 1e-8 * l1

        # do SGD step
        loss.backward()
        # update gradients on all the GPUs
        average_gradients(model, args)
        optimizer.step()

        # measure elapsed time
        batch_time.update(time.time() - end)
        end = time.time()

        if i % args.print_freq == 0 or i == len(train_loader) - 1:
            logger.info(
                "Epoch: %d, Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[ACC:%.3f, rnaSEQ:%.3f, TFs:%.3f], roc[ACC: %.3f, rnSEQ:%.3f, TFs: %.3f]" % (
                    epoch + 1, i, len(train_loader) - 1,
                    loss, perc_uncertainty, custom_accuracy[0], custom_accuracy[1], custom_accuracy[2], aucs[0], aucs[1], aucs[2]))


def validate(val_loader, model, criterion, args, logger):
    batch_time = AverageMeter('Time', ':6.3f')
    losses = AverageMeter('Loss', ':.4e')

    # switch to evaluate mode
    total_ACC, total_RNA, total_TFs = 0, 0, 0
    model.eval()
    
    perc_uncertainty = 0.0

    with torch.no_grad():
        end = time.time()
        for i, (input, target) in enumerate(val_loader):
            if args.gpu is not None:
                input = input.cuda(args.gpu, non_blocking=True)
                target = target.cuda(args.gpu, non_blocking=True)

            # compute output
            # output = model(input, args)
            output = model(input)
            loss = criterion(output, target)

            # measure accuracy and record loss
            # acc = getAccuracy(output, target, args)
            p = find_perc_uncertainty(output, 0.4, 0.6)
            custom_accuracy = getCustomAccuracy2(output, target, args)
            aucs = getAUCscore(output, target, args, logger)

            total_ACC += np.median(aucs[0])
            total_RNA += np.median(aucs[1])
            total_TFs += np.median(aucs[2])

            losses.update(loss.item(), input.size(0))
            # top1.update(acc, input.size(0))

            # measure elapsed time
            batch_time.update(time.time() - end)
            end = time.time()

            if i % args.print_freq == 0 or i == len(val_loader) - 1:
                # logger.info("batch: %d, loss: %.3f; valid accuracy: ACC marks: %.3f, RNA-seq: %.3f, TFs: %.3f" % (i+1, loss, aucs[0], aucs[1], aucs[2]))
                logger.info(
                    "Batch: %d, Loss: %.3f, custom[ACC: %.3f, rnaSEQ: %.3f, TFs: %.3f], roc[ACC: %.3f, rnSEQ: %.3f, TFs: %.3f]" % (
                        i, loss, custom_accuracy[0], custom_accuracy[1], custom_accuracy[2], aucs[0], aucs[1], aucs[2]))
            perc_uncertainty += p

                # logger.info("batch: %d, valid accuracy: %.3f" % (i+1, acc))
        logger.info("percentage of uncertainty in training prediction: {}".format(perc_uncertainty / len(val_loader)))


    total_ACC /= len(val_loader)
    total_RNA /= len(val_loader)
    total_TFs /= len(val_loader)

    return np.median([total_ACC, total_RNA, total_TFs])


def average_gradients(model, args):
    """ Gradient averaging. """
    # size = float(dist.get_world_size())
    size = float(args.world_size)
    for param in model.parameters():
        dist.all_reduce(param.grad.data, op=dist.ReduceOp.SUM)
        param.grad.data /= size

        
def save_checkpoint(state, is_best, modelName):
    torch.save(state, modelName + "_checkpoint.pth.tar")
    if is_best:
        shutil.copyfile(modelName + "_checkpoint.pth.tar", modelName + '_model_best.pth.tar')


class AverageMeter(object):
    """Computes and stores the average and current value"""
    def __init__(self, name, fmt=':f'):
        self.name = name
        self.fmt = fmt
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count

    def __str__(self):
        fmtstr = '{name} {val' + self.fmt + '} ({avg' + self.fmt + '})'
        return fmtstr.format(**self.__dict__)


def adjust_learning_rate(optimizer, epoch, args, lr_scheduler):
    """Sets the learning rate to the initial LR decayed by 10 every 30 epochs"""
    # lr = args.w_lr * (0.1 ** (epoch // 30))   # option 1
    lr_scheduler.step()     # option 2
    lr = lr_scheduler.get_lr()[0]

    for param_group in optimizer.param_groups:
        param_group['lr'] = lr

if __name__ == '__main__':
    main()
