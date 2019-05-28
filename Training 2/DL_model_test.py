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
import torch.tensor

# Setting arguments
parser = argparse.ArgumentParser("deepBrain2")
parser.add_argument('--name', required=True)
parser.add_argument('--DataDir', required=True)
parser.add_argument('--TrainingDataFile', required=True)
parser.add_argument('--TrainingLabelFile', required=True)
parser.add_argument('--TestingDataFile', required=True)
parser.add_argument('--TestingLabelFile', required=True)
parser.add_argument('--w_lr', type=float, default=1e-2, help='lr for weights')
parser.add_argument('--w_lr_min', type=float, default=8e-7, help='minimum lr for weights')
parser.add_argument('--w_momentum', type=float, default=0.9)
parser.add_argument('--w_weight_decay', type=float, default=5e-7, help='L2 Regularization Coefficient')
parser.add_argument('--w_grad_clip', type=float, default=5.,
                    help='gradient clipping for weights')
parser.add_argument('--print_freq', type=int, default=1, help='print frequency')
parser.add_argument('--init_channels', type=int, default=16)
parser.add_argument('--layers', type=int, default=8)
parser.add_argument('--BATCH_SIZE', type=int, default=16, metavar='N',
                    help='mini-batch size, this is the total batch size of all GPUs on the current '
                         'node when using Data Parallel or Distributed Data Parallel')
parser.add_argument('--seed', type=int, default=0, metavar='S',
                        help='random seed (default: 1)')
parser.add_argument('--workers', type=int, default=4, help='# of workers')
parser.add_argument('--alpha_lr', type=float, default=3e-4, help='lr for alpha')
parser.add_argument('--alpha_weight_decay', type=float, default=1e-4, help='weight decay for alpha')
parser.add_argument('--world-size', default=-1, type=int, help='number of nodes for distributed training')
parser.add_argument('--rank', default=0, type=int, help='node rank for distributed training')
parser.add_argument('--dist-url', default='env://', type=str, help='url used to set up distributed training')
parser.add_argument('--dist-backend', default='nccl', type=str, help='distributed backend')
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

# another DeepSEA implementation 
# (copied from "https://github.com/kipoi/models/blob/master/DeepSEA/model_architecture.py")
# Note: didn't perform better than DeeperDeepSEA on the small data (tempTrain5.dat/chr5)
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




def getCustomAccuracy(predicted, target, args):
    # predicted = torch.round(torch.sigmoid(predicted))
    # predicted = torch.round(predicted)
    n_digits = 3
    _predicted = torch.round(predicted * 10 ** n_digits) / (10 ** n_digits)
    __predicted = torch.round(_predicted)

    N = predicted.size(0) * args.NUM_OUTPUTS

    truePred = torch.sum(torch.eq(__predicted, target)).item()
    acc_val = truePred / N

    # print(torch.sum(torch.eq(target, torch.ones(target.shape))).item())
    # print(torch.sum(torch.eq(predicted, torch.ones(target.shape))).item())
    return acc_val


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


def getCustomAccuracy3(predicted, target, args):
    preds = []
    targets = []
    for i in range(10):
        o = F.log_softmax(torch.autograd.Variable(predicted), dim=1)
        t = torch.autograd.Variable(target)

        _, pred = torch.max(o, dim=1)
        preds.append(pred.data)
        targets.append(t.data)

    preds = torch.cat(preds)
    targets = torch.cat(targets)


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

    # print('Medican AUCs: Accetylation marks: %.3f, RNA-seq: %.3f, TFs: %.3f' % (np.median(aucs[:2]), np.median(aucs[2]), np.median(aucs[3:])))

    return np.median(aucs[:2]), np.median(aucs[2]), np.median(aucs[3:])


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
    output = output.detach().numpy()
    return np.sum(np.logical_and(output>=low, output<=high))/output.size    # returns the proportion of tensor elements are within the range


def train(train_loader, model, criterion, optimizer, epoch, args, logger):

    # switch to train mode
    model.train()
    # perc_uncertainty = 0.0

    for i, (input, target) in enumerate(train_loader):
        if (i == len(train_loader) - 1):
            x = 0
        # compute output
        # output = model(input, args)
        output = model(input)
        loss = criterion(output, target.squeeze(1))

        # measure accuracy and record loss
        # acc = getCustomAccuracy(output, target, args)
        perc_uncertainty = find_perc_uncertainty(output, 0.4, 0.6)
        custom_accuracy = getCustomAccuracy2(output, target, args)
        # tAccuracy = getCustomAccuracy3(output, target, args)
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
        optimizer.step()

        logger.info(
            "Epoch: %d, Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[ACC:%.3f, rnaSEQ:%.3f, TFs:%.3f], roc[ACC:%.3f, rnSEQ:%.3f, TFs:%.3f]" % (
            epoch + 1, i, len(train_loader) - 1,
            loss, perc_uncertainty, custom_accuracy[0], custom_accuracy[1], custom_accuracy[2], aucs[0], aucs[1], aucs[2]))



def validate(val_loader, model, criterion, args, logger):

    # switch to evaluate mode
    model.eval()
    total_ACC, total_RNA, total_TFs = 0, 0, 0
    perc_uncertainty = 0.0

    # losses = []

    with torch.no_grad():
        for i, (input, target) in enumerate(val_loader):
            # compute output
            # output = model(input, args)
            output = model(input)
            loss = criterion(output, target)
            # losses[i] = loss

            # measure accuracy and record loss
            # acc = getCustomAccuracy(output, target, args)
            p = find_perc_uncertainty(output, 0.4, 0.6)
            custom_accuracy = getCustomAccuracy2(output, target, args)
            aucs = getAUCscore(output, target, args, logger)

            total_ACC += np.median(aucs[0])
            total_RNA += np.median(aucs[1])
            total_TFs += np.median(aucs[2])

            if i % args.print_freq == 0 or i == len(val_loader) - 1:
                # progress._print(i)
                # logger.info("batch: %d, loss: %.3f; valid accuracy: custom_accuracy_metric: %.3f, ACC marks: %.3f, RNA-seq: %.3f, TFs: %.3f" % (i+1, loss, acc, aucs[0], aucs[1], aucs[2]))
                logger.info("Batch: %d/%d, Loss: %.3f, perc_uncertainty: %.3f, custom[ACC:%.3f, rnaSEQ:%.3f, TFs:%.3f], roc[ACC: %.3f, rnSEQ: %.3f, TFs:%.3f]" % (
                        i, len(val_loader)-1, loss, perc_uncertainty, custom_accuracy[0], custom_accuracy[1], custom_accuracy[2], aucs[0], aucs[1], aucs[2]))
            perc_uncertainty += p

        # logger.info(' * Acc@1 {top1.avg:.3f}'.format(top1=acc))
        logger.info("percentage of uncertainty in validation prediction: {}".format(perc_uncertainty / len(val_loader)))

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
    args = parser.parse_args()
    use_cuda = torch.cuda.is_available()
    torch.manual_seed(args.seed)
    device = torch.device("cuda" if use_cuda else "cpu")

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
    valLoader = tdata.DataLoader(dataset=valDataset, batch_size=args.BATCH_SIZE, shuffle=False,
                                 num_workers=args.workers, pin_memory=True)

    # build the model and criterion
    # model = BuildModel(args).to(device)
    # model = DeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS).to(device)
    model = DeeperDeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS).to(device)
    # model = get_model(load_weights=False)

    # Add a sigmoid activation function to the output.  Use a binary cross entropy
    # criterion = nn.BCEWithLogitsLoss()
    criterion = nn.CrossEntropyLoss()

    # optimiser = topti.Adam(model.parameters(), lr=args.w_lr)  # Minimise the loss using the Adam algorithm.
    # optimiser = torch.optim.Adam(model.parameters(), args.w_lr, betas=(0.5, 0.999),
    #                                weight_decay=args.w_weight_decay)
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
        train(trainLoader, model, criterion, optimiser, epoch, args, logger)

        # evaluate on validation set
        acc1 = validate(valLoader, model, criterion, args, logger)

        # scheduler.step(acc1[1])

        # remember best acc@1 and save checkpoint
        is_best = acc1 > best_acc1
        best_acc1 = max(acc1, best_acc1)

    # if (args.save_model):
    #     torch.save(model.state_dict(), "mnist_cnn.pt")
    torch.save(model.state_dict(), "deepBrain2_cnn.pt")


if __name__ == '__main__':
    main()
