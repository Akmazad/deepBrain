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
import utils
from models.search_cnn import SearchCNN
from architect import Architect
import sys


parser = argparse.ArgumentParser("deepBrain2")
parser.add_argument('--name', required=True)
parser.add_argument('--DataDir', required=True)
parser.add_argument('--TrainingDataFile', required=True)
parser.add_argument('--TrainingLabelFile', required=True)
parser.add_argument('--TestingDataFile', required=True)
parser.add_argument('--TestingLabelFile', required=True)
parser.add_argument('--w_lr', type=float, default=0.01, help='lr for weights')
parser.add_argument('--w_lr_min', type=float, default=0.001, help='minimum lr for weights')
parser.add_argument('--w_momentum', type=float, default=0.9)
parser.add_argument('--w_weight_decay', type=float, default=3e-4)
parser.add_argument('--w_grad_clip', type=float, default=5.,
                    help='gradient clipping for weights')
parser.add_argument('--print_freq', type=int, default=1, help='print frequency')
parser.add_argument('--init_channels', type=int, default=16)
parser.add_argument('--layers', type=int, default=8)
parser.add_argument('--BATCH_SIZE', type=int, default=256, metavar='N',
                    help='mini-batch size, this is the total batch size of all GPUs on the current '
                         'node when using Data Parallel or Distributed Data Parallel')
parser.add_argument('--seed', type=int, default=1, metavar='S',
                        help='random seed (default: 1)')
parser.add_argument('--workers', type=int, default=1, help='# of workers')
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
parser.add_argument('--KERNEL_SIZE', type=int, default=4)
parser.add_argument('--POOLING_TH', type=int, default=2)  # max-pooling is used with 20%
parser.add_argument('--DROPOUT_l1', type=float, default=0.2)
parser.add_argument('--DROPOUT_l2', type=float, default=0.2)
parser.add_argument('--DROPOUT_l3', type=float, default=0.5)
parser.add_argument('--NUM_OUTPUTS', type=int, default=131)
parser.add_argument('--SEQ_LEN', type=int, default=1000)

best_acc1 = 0


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


def getAccuracy(predicted, target, args):
    # predicted = torch.round(torch.sigmoid(predicted))
    predicted = torch.round(predicted)
    N = predicted.size(0) * args.NUM_OUTPUTS

    truePred = torch.sum(torch.eq(predicted, target)).item()
    acc_val = truePred / N

    print(torch.sum(torch.eq(target, torch.ones(target.shape))).item())
    print(torch.sum(torch.eq(predicted, torch.ones(target.shape))).item())
    return acc_val


def getAUCscore(predicted, target, args, logger):
    # _predicted = torch.round(predicted).detach()
    _predicted = predicted.detach()
    _target = target.detach()

    aucs = np.zeros(args.NUM_OUTPUTS, dtype=np.float)
    for i in range(args.NUM_OUTPUTS):
        try:
            auc = roc_auc_score(_target.cpu().numpy()[:, i], _predicted.cpu().numpy()[:, i])
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


def train(train_loader, valid_loader, model, architect, w_optim, alpha_optim, lr, epoch, logger, args):
    stats = utils.SumMeter()
    # logger.info('8.1')
    # top1 = utils.AverageMeter()
    # top5 = utils.AverageMeter()
    losses = utils.AverageMeter()
    # logger.info('8.2')

    cur_step = epoch * len(train_loader)

    # logger.info('8.3')

    model.train()
    # logger.info('8.4')

    train_iter = iter(train_loader)
    valid_iter = iter(valid_loader)
    # logger.info('8.5')

    for step, ((input, target), (val_X, val_y)) in enumerate(zip(train_loader, valid_loader)):
        # if not torch.cuda.is_available():
        #     logging.info('no gpu device available')
        #     sys.exit(1)

        # if args.gpu is not None:
        #     input = input.cuda(args.gpu, non_blocking=True)
        #     val_X = val_X.cuda(args.gpu, non_blocking=True)
        # target = target.cuda(args.gpu, non_blocking=True)
        # val_y = val_y.cuda(args.gpu, non_blocking=True)

        logger.info('8.8')
        # logger.info('8.9')

        # phase 2. architect step (alpha)
        alpha_optim.zero_grad()
        logger.info('8.10')
        try:
            architect.unrolled_backward(input, target, val_X, val_y, lr, w_optim, logger)
        except Exception as e:
            logger.info("problem: {}".format(e))
        logger.info('8.11')
        alpha_optim.step()
        logger.info('8.12')

        # phase 1. child network step (w)
        w_optim.zero_grad()
        # logger.info('8.13')
        output = model(input)
        logger.info('8.14')
        loss = model.criterion(output, target)
        logger.info('8.15')
        loss.backward()
        logger.info('8.16')
        # gradient clipping
        nn.utils.clip_grad_norm_(model.weights(), args.w_grad_clip)
        # logger.info('8.17')
        w_optim.step()
        # logger.info('8.18')

        acc = getAccuracy(output, target, args)
        # logger.info('8.19')

        # logger.info('8.20')

        # logger.info('8.21')

        if step % args.print_freq == 0 or step == len(train_loader) - 1:
            # progress._print_(i)
            logger.info("Epoch: %d, Batch: %d/%d, Loss: %.3f, acc: %.3f" % (epoch + 1, step, len(train_loader) - 1,
                                                                            loss, acc))

        cur_step += 1

    # logger.info("Train: [{:2d}/{}] Final Prec@1 {:.4}".format(epoch + 1, args.epochs, stats.MCC()))



def validate(valid_loader, model, epoch, cur_step, logger, args):
    top1 = AverageMeter('Acc@1', ':6.2f')
    # set validation mode
    model.eval()

    with torch.no_grad():
        for step, (input, target) in enumerate(valid_loader):
            if args.gpu is not None:
                input = input.cuda(args.gpu, non_blocking=True)
            target = target.cuda(args.gpu, non_blocking=True)

            output = model(input)
            loss = model.criterion(output, target)

            acc = getAccuracy(output, target, args)
            top1.update(acc, input.size(0))

            if step % args.print_freq == 0 or step == len(valid_loader) - 1:
                # progress._print(i)
                logger.info("batch: %d, losses: %.3f,  valid accuracy: %.3f" % (step + 1, loss, acc))

        logger.info(' * Acc@1 {top1.avg:.3f}'.format(top1=top1))

    return top1.avg


def adjust_learning_rate(optimizer, epoch, args):
    """Sets the learning rate to the initial LR decayed by 10 every 30 epochs"""
    lr = args.w_lr * (0.1 ** (epoch // 30))
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr


def main():
    # Training settings
    args = parser.parse_args()
    use_cuda = torch.cuda.is_available()
    torch.manual_seed(args.seed)
    device = torch.device("cuda" if use_cuda else "cpu")

    DataPath = args.DataDir
    logger = get_logger(os.path.join(os.path.join(DataPath, args.name, 'DL_model_test_DARTS'), "{}.log".format(args.name)))

    net_crit = nn.BCELoss()
    model = SearchCNN(args.CONV1_INPUT_CHANNELS, args.init_channels, args.NUM_OUTPUTS, args.layers,
                      net_crit)  # DARTS model

    # weights optimizer
    w_optim = torch.optim.SGD(model.weights(), args.w_lr, momentum=args.w_momentum,
                              weight_decay=args.w_weight_decay)
    # alphas optimizer
    alpha_optim = torch.optim.Adam(model.alphas(), args.alpha_lr, betas=(0.5, 0.999),
                                   weight_decay=args.alpha_weight_decay)

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
    # model = DeeperDeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS).to(device)


    # Add a sigmoid activation function to the output.  Use a binary cross entropy
    # criterion = nn.BCEWithLogitsLoss()
    # criterion = nn.BCELoss()

    # optimiser = topti.Adam(model.parameters(), lr=args.w_lr)  # Minimise the loss using the Adam algorithm.

    # model = DeepSEA(args.SEQ_LEN, args.NUM_OUTPUTS).to(device)
    # criterion = nn.BCELoss()
    # optimiser = (torch.optim.SGD, {"lr": args.w_lr, "weight_decay": 1e-6, "momentum": 0.9})

    lr_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        w_optim, args.nEpochs, eta_min=args.w_lr_min)
    # logger.info('6')

    architect = Architect(model, args.w_momentum, args.w_weight_decay)

    best_acc1 = 0.0
    for epoch in range(args.start_epoch, args.nEpochs):

        # adjust_learning_rate(optimiser, epoch, args)
        lr_scheduler.step()
        lr = lr_scheduler.get_lr()[0]

        train(trainLoader, valLoader, model, architect, w_optim, alpha_optim, lr, epoch, logger, args)
        logger.info('training done')

        # evaluate on validation set
        # acc1 = validate(valLoader, model, criterion, args, logger)
        cur_step = (epoch + 1) * len(trainLoader)
        top1 = validate(valLoader, model, epoch, cur_step, logger, args)
        logger.info('validation done')
        # genotype
        genotype = model.genotype()
        logger.info("genotype = {}".format(genotype))

        # remember best acc@1 and save checkpoint
        is_best = top1 > best_acc1
        best_acc1 = max(top1, best_acc1)

    # if (args.save_model):
    #     torch.save(model.state_dict(), "mnist_cnn.pt")
    torch.save(model.state_dict(), "deepBrain2_cnn.pt")

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

if __name__ == '__main__':
    main()