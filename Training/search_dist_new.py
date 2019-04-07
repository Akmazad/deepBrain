import argparse
import os
import random
import shutil
import time
import warnings
import sys

import torch
# import torch.nn as nn
import torch.nn.parallel
import torch.backends.cudnn as cudnn
# import torch.distributed as dist
import torch.optim
import torch.multiprocessing as mp
import torch.utils.data
import torch.utils.data.distributed
import torchvision.transforms as transforms
import torchvision.datasets as datasets
import torchvision.models as models

""" Search cell """
import os
import torch
import torch.nn as nn
import numpy as np
from tensorboardX import SummaryWriter
from config import SearchConfig
import utils
from models.search_cnn import SearchCNN
from architect import Architect
from visualize import plot

from random import shuffle
import pickle

# special import starts: this imports are helpful (may be) for multi-node multi-GPU support
from torch.cuda.comm import broadcast_coalesced
from torch.cuda import nccl
import torch.distributed as dist

if dist.is_available():
    from torch.distributed.distributed_c10d import _get_default_group

# from ..modules import Module
# from .replicate import replicate
# from .scatter_gather import scatter_kwargs, gather
# from .parallel_apply import parallel_apply
from torch.cuda._utils import _get_device_index

# special import ends


config = SearchConfig()

logger = utils.get_logger(os.path.join(config.path, "{}.log".format(config.name)))


# config.print_params(logger.info)

def main():
    try:
        mp.set_start_method('spawn')
    except RuntimeError:
        pass

    logger.info("Logger is set - training start")

    if config.seed is not None:
        random.seed(config.seed)
        torch.manual_seed(config.seed)
        cudnn.deterministic = True
        warnings.warn('You have chosen to seed training. '
                      'This will turn on the CUDNN deterministic setting, '
                      'which can slow down your training considerably! '
                      'You may see unexpected behavior when restarting '
                      'from checkpoints.')
    logger.info('1')

    if config.gpu is not None:
        warnings.warn('You have chosen a specific GPU. This will completely '
                      'disable data parallelism.')

    if config.dist_url == "env://" and config.world_size == -1:
        config.world_size = int(os.environ["WORLD_SIZE"])

    logger.info("config.world_size: {}".format(config.world_size))

    config.distributed = config.world_size > 1 or config.multiprocessing_distributed
    logger.info("config.distributed: {}".format(config.distributed))

    ngpus_per_node = torch.cuda.device_count()
    if config.multiprocessing_distributed:
        # Since we have ngpus_per_node processes per node, the total world_size
        # needs to be adjusted accordingly
        config.world_size = ngpus_per_node * config.world_size
        # Use torch.multiprocessing.spawn to launch distributed processes: the
        # main_worker process function
        mp.spawn(main_worker, nprocs=ngpus_per_node, args=(ngpus_per_node, config))
    else:
        # Simply call main_worker function
        main_worker(config.gpu, ngpus_per_node, config)


def main_worker(gpu, ngpus_per_node, args):
    global best_acc1
    args.gpu = gpu
    logger.info("args.gpu: {}".format(args.gpu))

    if args.gpu is not None:
        logger.info("Use GPU: {} for training".format(args.gpu))

    if args.distributed:
        if args.dist_url == "env://" and args.rank == -1:
            args.rank = int(os.environ["RANK"])
        if args.multiprocessing_distributed:
            # For multiprocessing distributed training, rank needs to be the
            # global rank among all the processes
            args.rank = args.rank * ngpus_per_node + args.gpu
        dist.init_process_group(backend=args.dist_backend, init_method=args.dist_url,
                                world_size=args.world_size, rank=args.rank)

    logger.info('2')
    # get data with meta info
    input_size, input_channels, n_classes, train_data = utils.get_data(
        config.train_data, config.train_label, config.data_path, logger, cutout_length=0, validation=False)

    # net_crit = nn.BCEWithLogitsLoss().to(device)
    net_crit = nn.BCEWithLogitsLoss().cuda(args.gpu)
    # create model
    model = SearchCNN(input_channels, config.init_channels, n_classes, config.layers, net_crit)

    logger.info('3')
    if args.distributed:
        # For multiprocessing distributed, DistributedDataParallel constructor
        # should always set the single device scope, otherwise,
        # DistributedDataParallel will use all available devices.
        if args.gpu is not None:
            logger.info('3.1.0')
            torch.cuda.set_device(args.gpu)
            logger.info('3.1.1')
            model.cuda(args.gpu)
            logger.info('3.1.2')
            # When using a single GPU per process and per
            # DistributedDataParallel, we need to divide the batch size
            # ourselves based on the total number of GPUs we have
            args.batch_size = int(args.batch_size / ngpus_per_node)
            logger.info('3.1.3')
            args.workers = int(args.workers / ngpus_per_node)
            logger.info('3.1.4')
            model = torch.nn.parallel.DistributedDataParallel(model, device_ids=[args.gpu], output_device=args.gpu)
            logger.info('3.1.5')
        else:
            logger.info('3.2.0')
            model.cuda()
            logger.info('3.2.1')
            # DistributedDataParallel will divide and allocate batch_size to all
            # available GPUs if device_ids are not set
            model = torch.nn.parallel.DistributedDataParallel(model)
            logger.info('3.2.2')
    elif args.gpu is not None:
        logger.info('3.3.0')
        torch.cuda.set_device(args.gpu)
        logger.info('3.3.1')
        model = model.cuda(args.gpu)
        logger.info('3.3.2')
    else:
        # DataParallel will divide and allocate batch_size to all available GPUs
        logger.info('3.4.0')
        model = torch.nn.DataParallel(model).cuda()
        logger.info('3.4.1')
    logger.info('4')

    # weights optimizer
    w_optim = torch.optim.SGD(model.module.weights(), config.w_lr, momentum=config.w_momentum,
                              weight_decay=config.w_weight_decay)
    alpha_optim = torch.optim.Adam(model.module.alphas(), config.alpha_lr, betas=(0.5, 0.999),
                                   weight_decay=config.alpha_weight_decay)

    # optionally resume from a checkpoint
    if args.resume:
        if os.path.isfile(args.resume):
            logger.info("=> loading checkpoint '{}'".format(args.resume))
            checkpoint = torch.load(args.resume)
            args.start_epoch = checkpoint['epoch']
            best_acc1 = checkpoint['best_acc1']
            if args.gpu is not None:
                # best_acc1 may be from a checkpoint from a different GPU
                best_acc1 = best_acc1.to(args.gpu)
            model.load_state_dict(checkpoint['state_dict'])
            w_optim.load_state_dict(checkpoint['optimizer'])
            logger.info("=> loaded checkpoint '{}' (epoch {})".format(args.resume, checkpoint['epoch']))
        else:
            logger.info("=> no checkpoint found at '{}'".format(args.resume))

    cudnn.benchmark = True

    logger.info('5')

    # # Data loading code
    # n_train = len(train_data)
    # split = int(0.8 * n_train)
    # indices = list(range(n_train))
    # shuffle(indices)
    # trainIndices = indices[:split]
    # testIndices = indices[split:]
    # with open(config.data_path + "trainTestIndices.pickle", "wb") as indicesFile:
    #     pickle.dump(trainIndices, indicesFile)
    #     pickle.dump(testIndices, indicesFile)
    #
    # with open(config.data_path + "trainTestIndices.pickle", "rb") as indicesFile:
    #     trainIndices = pickle.load(indicesFile)
    # n_train = len(trainIndices)
    # split = n_train // 2
    # train_sampler = torch.utils.data.sampler.SubsetRandomSampler(trainIndices[:split])
    # valid_sampler = torch.utils.data.sampler.SubsetRandomSampler(trainIndices[split:])

    if args.distributed:
        train_sampler = torch.utils.data.distributed.DistributedSampler(train_data)
        # this is for checking purpose, must be changed
        valid_sampler = torch.utils.data.distributed.DistributedSampler(train_data)
    else:
        train_sampler = None
        valid_sampler = None

    train_loader = torch.utils.data.DataLoader(train_data,
                                               batch_size=config.batch_size,
                                               sampler=train_sampler,
                                               num_workers=config.workers,
                                               pin_memory=False)
    valid_loader = torch.utils.data.DataLoader(train_data,
                                               batch_size=config.batch_size,
                                               sampler=valid_sampler,
                                               num_workers=config.workers,
                                               pin_memory=False)

    logger.info("current rank: ".format(args.rank))

    # if args.evaluate:
    #     validate(valid_loader, model, net_crit, args)
    #     return
    lr_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(w_optim, config.epochs, eta_min=config.w_lr_min)
    logger.info('6')
    architect = Architect(model, config.w_momentum, config.w_weight_decay)

    # training loop
    best_genotype = None
    best_top1 = 0
    for epoch in range(args.start_epoch, args.epochs):
        # if args.distributed:
        #     train_sampler.set_epoch(epoch)

        lr_scheduler.step()
        lr = lr_scheduler.get_lr()[0]

        logger.info('7')

        model.module.print_alphas()
        logger.info('8')

        # training
        train(train_loader, valid_loader, model, architect, w_optim, alpha_optim, lr, epoch, args)
        logger.info('9')

        # validation
        cur_step = (epoch + 1) * len(train_loader)
        top1 = validate(valid_loader, model, net_crit, args)
        logger.info('10')

        # log
        # genotype
        genotype = model.module.genotype()
        logger.info("genotype = {}".format(genotype))

        # genotype as a image
        plot_path = os.path.join(config.plot_path, "EP{:02d}".format(epoch + 1))
        caption = "Epoch {}".format(epoch + 1)
        plot(genotype.normal, plot_path + "-normal", caption)
        plot(genotype.reduce, plot_path + "-reduce", caption)

        # save
        is_best = None
        if best_top1 < top1:
            best_top1 = top1
            best_genotype = genotype
            is_best = True
        else:
            is_best = False

        if not args.multiprocessing_distributed or (args.multiprocessing_distributed
                                                    and args.rank % ngpus_per_node == 0):
            save_checkpoint({
                'epoch': epoch + 1,
                'arch': architect,
                'state_dict': model.state_dict(),
                'best_acc1': best_top1,
                'optimizer': w_optim.state_dict(),
            }, is_best, config.path)

    logger.info("Final best Prec@1 = {:.4%}".format(best_top1))
    logger.info("Best Genotype = {}".format(best_genotype))


def train(train_loader, valid_loader, model, architect, w_optim, alpha_optim, lr, epoch, args):
    batch_time = AverageMeter()
    data_time = AverageMeter()
    losses = AverageMeter()
    top1 = AverageMeter()
    top5 = AverageMeter()

    cur_step = epoch * len(train_loader)

    # switch to train mode
    model.train()

    end = time.time()

    train_iter = iter(train_loader)
    valid_iter = iter(valid_loader)
    logger.info('8.1')

    for step, ((trn_X, trn_y), (val_X, val_y)) in enumerate(zip(train_loader, valid_loader)):
        # measure data loading time
        data_time.update(time.time() - end)
        # pass a mini-batch of data into GPUs if specified
        if args.gpu is not None:
            trn_X, trn_y = trn_X.cuda(args.gpu, non_blocking=True), trn_y.cuda(args.gpu, non_blocking=True)
            val_X, val_y = val_X.cuda(args.gpu, non_blocking=True), val_y.cuda(args.gpu, non_blocking=True)
        N = trn_X.size(0)

        # Architect step (alpha)
        alpha_optim.zero_grad()
        architect.unrolled_backward(trn_X, trn_y, val_X, val_y, lr, w_optim, logger)
        alpha_optim.step()

        # Child network step (w)
        w_optim.zero_grad()
        output = model(trn_X)  # compute output
        loss = model.module.criterion(output, trn_y)
        loss.backward()
        nn.utils.clip_grad_norm_(model.module.weights(), config.w_grad_clip)  # gradient clipping
        w_optim.step()

        # measure accuracy and record loss
        acc1, acc5 = accuracy(output, trn_y, topk=(1, 5))
        losses.update(loss.item(), trn_X.size(0))
        top1.update(acc1[0], trn_X.size(0))
        top5.update(acc5[0], trn_X.size(0))

        # # compute gradient and do SGD step
        # optimizer.zero_grad()
        # loss.backward()
        # optimizer.step()
        #
        logger.info('8.2')

        if step % config.print_freq == 0 or step == len(train_loader) - 1:
            logger.info(
                "Train: [{:2d}/{}] Step {:03d}/{:03d} Loss {losses.avg:.3f} "
                "Prec@(1,5) ({top1:.1%}, {top5:.3})".format(
                    epoch + 1, config.epochs, step, len(train_loader) - 1, losses=losses,
                    top1=top1, top5=top5))

        cur_step += 1

    logger.info("Train: [{:2d}/{}] Final Prec@1 {:.4}".format(epoch + 1, config.epochs, top5))


def validate(val_loader, model, criterion, args):
    batch_time = AverageMeter()
    losses = AverageMeter()
    top1 = AverageMeter()
    top5 = AverageMeter()

    # switch to evaluate mode
    model.eval()

    with torch.no_grad():
        end = time.time()
        for i, (input, target) in enumerate(val_loader):
            if args.gpu is not None:
                input = input.cuda(args.gpu, non_blocking=True)
            target = target.cuda(args.gpu, non_blocking=True)

            # compute output
            output = model(input)
            loss = criterion(output, target)

            # measure accuracy and record loss
            acc1, acc5 = accuracy(output, target, topk=(1, 5))
            losses.update(loss.item(), input.size(0))
            top1.update(acc1[0], input.size(0))
            top5.update(acc5[0], input.size(0))

            # measure elapsed time
            batch_time.update(time.time() - end)
            end = time.time()

            if i % args.print_freq == 0:
                logger.info('Test: [{0}/{1}]\t'
                            'Time {batch_time.val:.3f} ({batch_time.avg:.3f})\t'
                            'Loss {loss.val:.4f} ({loss.avg:.4f})\t'
                            'Acc@1 {top1.val:.3f} ({top1.avg:.3f})\t'
                            'Acc@5 {top5.val:.3f} ({top5.avg:.3f})'.format(
                    i, len(val_loader), batch_time=batch_time, loss=losses,
                    top1=top1, top5=top5))

        logger.info(' * Acc@1 {top1.avg:.3f} Acc@5 {top5.avg:.3f}'
                    .format(top1=top1, top5=top5))

    return top1.avg


def save_checkpoint(state, is_best, ckpt_dir, filename='checkpoint.pth.tar'):
    filename = os.path.join(ckpt_dir, 'checkpoint.pth.tar')
    torch.save(state, filename)
    if is_best:
        best_filename = os.path.join(ckpt_dir, 'best.pth.tar')
        shutil.copyfile(filename, best_filename)


class AverageMeter(object):
    """Computes and stores the average and current value"""

    def __init__(self):
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


def adjust_learning_rate(optimizer, epoch, args):
    """Sets the learning rate to the initial LR decayed by 10 every 30 epochs"""
    lr = args.lr * (0.1 ** (epoch // 30))
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr


def accuracy(output, target, topk=(1,)):
    """Computes the accuracy over the k top predictions for the specified values of k"""
    with torch.no_grad():
        maxk = max(topk)
        batch_size = target.size(0)

        _, pred = output.topk(maxk, 1, True, True)
        pred = pred.t()
        correct = pred.eq(target.view(1, -1).expand_as(pred))

        res = []
        for k in topk:
            correct_k = correct[:k].view(-1).float().sum(0, keepdim=True)
            res.append(correct_k.mul_(100.0 / batch_size))
        return res


if __name__ == '__main__':
    main()
