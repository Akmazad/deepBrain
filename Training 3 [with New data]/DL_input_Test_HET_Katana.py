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

import torch.tensor

class config():
  def __init__(self):
    self.name = "deepbrainStaticConvnet"
    self.DataDir = "/srv/scratch/z3526914/DeepBrain/Data/"
#     self.TrainingDataFile = "tempTrain5_trainingData_TF_filtered_chr5_value.npy"
#     self.TrainingLabelFile = "tempTrain5_trainingData_TF_filtered_chr5_label_all.npy"
#     self.TestingDataFile = "tempTrain5_validationData_TF_filtered_chr5_value.npy"
#     self.TestingLabelFile = "tempTrain5_validationData_TF_filtered_chr5_label_all.npy"

    self.TrainingDataFile = "HumanFC_ENCODE_EpiMap_nonZero.bin.Seq_Labels_trainingData_chr1_value.npy"
    self.TrainingLabelFile = "HumanFC_ENCODE_EpiMap_nonZero.bin.Seq_Labels_trainingData_chr1_label_all.npy"
    self.TestingDataFile = "HumanFC_ENCODE_EpiMap_nonZero.bin.Seq_Labels_validationData_chr1_value.npy"
    self.TestingLabelFile = "HumanFC_ENCODE_EpiMap_nonZero.bin.Seq_Labels_validationData_chr1_label_all.npy"
    
    self.w_lr = 1e-2
    self.w_lr_min = 8e-7
    self.w_momentum = 0.9
    self.w_weight_decay = 5e-7
    self.w_grad_clip = 5
    self.print_freq = 50
    self.init_channels = 16
    self.layers = 8
    self.BATCH_SIZE = 128
    self.val_batch_size = 128
    self.seed = 0
    self.workers = 4 
    self.alpha_lr = 3e-4
    self.alpha_weight_decay = 1e-4
    self.world_size = -1
    self.rank = 0
    self.dist_url = 'env://'
    self.dist_backend = 'nccl'
    self.gpu = None
    self.multiprocessing_distributed = True
    self.nEpochs = 15
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

# model zone
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
        # x = tfunc.batch_norm(x, torch.zeros(args.CONV1_OUTPUT_CHANNELS), torch.ones(args.CONV1_OUTPUT_CHANNELS))
        x = tfunc.dropout(x, args.DROPOUT_l1)

        x = self.conv_layer2(x)
        x = tfunc.leaky_relu(x, 0.01)
        x = tfunc.max_pool1d(x, args.POOLING_TH)
        # x = tfunc.batch_norm(x, torch.zeros(args.CONV2_OUTPUT_CHANNELS), torch.ones(args.CONV2_OUTPUT_CHANNELS))
        x = tfunc.dropout(x, args.DROPOUT_l2)

        x = self.conv_layer3(x)
        x = tfunc.leaky_relu(x, 0.1)
        x = tfunc.max_pool1d(x, args.POOLING_TH)
        # x = tfunc.batch_norm(x, torch.zeros(args.CONV3_OUTPUT_CHANNELS), torch.ones(args.CONV3_OUTPUT_CHANNELS))
        x = tfunc.dropout(x, args.DROPOUT_l3)

        # for the fully connected layer
        x = x.view(x.shape[0], -1)  # Flatten tensor.
        x = self.fc1(x)
        x = tfunc.leaky_relu(x, 0.01)
        #         x = tfunc.dropout(x, 0.2)
        # x = tfunc.batch_norm(x, torch.zeros(K).cuda(), torch.ones(K).cuda())
        x = torch.sigmoid(x)
        # x = F.log_softmax(x)

        return x

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
            # nn.ReLU(inplace=True),
            nn.Threshold(0, 1e-06),

            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.Dropout(p=0.2),

            nn.Conv1d(320, 480, kernel_size=conv_kernel_size),
            # nn.ReLU(inplace=True),
            nn.Threshold(0, 1e-06),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.Dropout(p=0.2),

            nn.Conv1d(480, 960, kernel_size=conv_kernel_size),
            # nn.ReLU(inplace=True),
            nn.Threshold(0, 1e-06),
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
            # nn.ReLU(inplace=True),
            nn.Threshold(0, 1e-06),
            nn.Linear(n_genomic_features, n_genomic_features),
            nn.Sigmoid())

    def forward(self, x):
        """Forward propagation of a batch.
        """
        out = self.conv_net(x)
        reshape_out = out.view(out.size(0), 960 * self.n_channels)
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

class ReCodeAlphabet(nn.Module):
    def __init__(self):
        super(ReCodeAlphabet, self).__init__()
        #
    def forward(self, input):
        # Swap ACGT to AGCT
        # array has shape (N, 4, 1, 1000)
        # pytorch doesn't support full indexing at the moment, at some point this should work: [:,:,torch.LongTensor([0,2,1,3])]
        input_reordered = [input[:,i,...] for i in [0,2,1,3]]
        input = torch.stack(input_reordered, dim=1)
        # slightly faster but ugly:
        #input = edit_tensor_in_numpy(input, lambda x: x[:,[0,2,1,3], ...])
        return input


class LambdaBase(nn.Sequential):
    def __init__(self, fn, *args):
        super(LambdaBase, self).__init__(*args)
        self.lambda_func = fn
        #
    def forward_prepare(self, input):
        output = []
        for module in self._modules.values():
            output.append(module(input))
        return output if output else input


class Lambda(LambdaBase):
    def forward(self, input):
        return self.lambda_func(self.forward_prepare(input))


def get_model(load_weights = True):
    deepsea_cpu = nn.Sequential( # Sequential,
        nn.Conv2d(4,320,(1, 8),(1, 1)),
        nn.Threshold(0, 1e-06),
        nn.MaxPool2d((1, 4),(1, 4)),
        nn.Dropout(0.2),
        nn.Conv2d(320,480,(1, 8),(1, 1)),
        nn.Threshold(0, 1e-06),
        nn.MaxPool2d((1, 4),(1, 4)),
        nn.Dropout(0.2),
        nn.Conv2d(480,960,(1, 8),(1, 1)),
        nn.Threshold(0, 1e-06),
        nn.Dropout(0.5),
        Lambda(lambda x: x.view(x.size(0),-1)), # Reshape,
        nn.Sequential(Lambda(lambda x: x.view(1,-1) if 1==len(x.size()) else x ),nn.Linear(50880,925)), # Linear,
        nn.Threshold(0, 1e-06),
        nn.Sequential(Lambda(lambda x: x.view(1,-1) if 1==len(x.size()) else x ),nn.Linear(925,919)), # Linear,
        nn.Sigmoid(),
    )
    if load_weights:
        deepsea_cpu.load_state_dict(torch.load('model_files/deepsea_cpu.pth'))
    return nn.Sequential(ReCodeAlphabet(), deepsea_cpu)


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
    n_digits = 3    # to have something like 0.499 = 0.5
    _predicted = torch.round(predicted * 10 ** n_digits) / (10 ** n_digits)
    __predicted = torch.round(_predicted)

    N = predicted.size(0)
    custom_accuracy = np.zeros(args.NUM_OUTPUTS, dtype=np.float)
    for i in range(args.NUM_OUTPUTS):
        truePred = torch.sum(torch.eq(__predicted[:, i], target[:, i])).item()
        custom_accuracy[i] = truePred/N

    return np.median(custom_accuracy[:288]), np.median(custom_accuracy[288:438]), np.median(custom_accuracy[438:])


def getAUCscore(predicted, target, args, logger):
    n_digits = 3
    _predicted = torch.round(predicted * 10**n_digits) / (10**n_digits)
    __predicted = torch.round(_predicted)

#     _predicted = torch.round(predicted).detach()
    __predicted = __predicted.detach()
    _target = target.detach()

    aucs = np.zeros(args.NUM_OUTPUTS, dtype=np.float)
    for i in range(args.NUM_OUTPUTS):
        try:
            auc = roc_auc_score(_target.cpu().numpy()[:, i], __predicted.cpu().numpy()[:, i])
            aucs[i] = auc
        except ValueError as e:
            pass
            # logger.info("NA (No positive (i.e. signal) in Test region)")

    return np.median(aucs[:288]), np.median(aucs[288:438]), np.median(aucs[438:])

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

        # measure accuracy and record loss
        # acc = getCustomAccuracy(output, target, args)
        perc_uncertainty = find_perc_uncertainty(output, 0.4, 0.6)
        custom_accuracy = getCustomAccuracy2(output, target, args)
        aucs = getAUCscore(output, target, args, logger)
        accs = getAccuracyScore(output, target, args, logger)
        f1s = getF1Score(output, target, args, logger)

        # compute gradient
        optimizer.zero_grad()

        # add l1 Sparsity
        l1 = 0
        for p in model.parameters():
            l1 = l1 + p.abs().sum()
        loss = loss + 1e-8 * l1

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
                "TRAINING: Epoch: %d, Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f], roc[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], accuracy[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], f1[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ]"
                % (epoch + 1, i, len(train_loader) - 1,
                   loss, perc_uncertainty, custom_accuracy[0],
                   custom_accuracy[1], custom_accuracy[2],
                   aucs[0], aucs[1], aucs[2],
                   accs[0], accs[1], accs[2],
                   f1s[0], f1s[1], f1s[2]))


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
            # losses[i] = loss

            # measure accuracy and record loss
            # acc = getCustomAccuracy(output, target, args)
            p = find_perc_uncertainty(output, 0.4, 0.6)
            custom_accuracy = getCustomAccuracy2(output, target, args)
            aucs = getAUCscore(output, target, args, logger)
            accs = getAccuracyScore(output, target, args, logger)
            f1s = getF1Score(output, target, args, logger)

            total_ACC += np.median(aucs[0])
            total_RNA += np.median(aucs[1])
            total_TFs += np.median(aucs[2])

            if i % args.print_freq == 0 or i == len(val_loader) - 1:
                # progress._print(i)
                # logger.info("batch: %d, loss: %.3f; valid accuracy: custom_accuracy_metric: %.3f, ACC marks: %.3f, RNA-seq: %.3f, TFs: %.3f" % (i+1, loss, acc, aucs[0], aucs[1], aucs[2]))
                #                 logger.info("VALIDATION: Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[HumanFC:%.3f, EpiMap:%.3f, TFs:%.3f], roc[HumanFC: %.3f, EpiMap: %.3f, TFs:%.3f]" % (i, len(val_loader)-1, loss, perc_uncertainty, custom_accuracy[0], custom_accuracy[1], custom_accuracy[2], aucs[0], aucs[1], aucs[2]))
                logger.info(
                    "VALIDATION: Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f], roc[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], accuracy[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ], f1[HumanFC: %.3f, EpiMap: %.3f, TFs: %.3f ]"
                    % (i, len(val_loader) - 1,
                       loss, p, custom_accuracy[0],
                       custom_accuracy[1], custom_accuracy[2],
                       aucs[0], aucs[1], aucs[2],
                       accs[0], accs[1], accs[2],
                       f1s[0], f1s[1], f1s[2]))

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
    logger = get_logger(os.path.join(os.path.join(DataPath, args.name, 'staticConvNet'), "{}.log".format(args.name)))

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
    model = DeeperDeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS).to(device)
    # model = get_model(load_weights=False)

    # Add a sigmoid activation function to the output.  Use a binary cross entropy
    criterion = nn.BCELoss().to(device)
    # criterion = nn.CrossEntropyLoss()

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
