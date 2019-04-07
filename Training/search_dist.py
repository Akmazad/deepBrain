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

device = torch.device("cuda:0")

# tensorboard
# writer = SummaryWriter(log_dir=os.path.join(config.path, "tb"))
# writer.add_text('config', config.as_markdown(), 0)

logger = utils.get_logger(os.path.join(config.path, "{}.log".format(config.name)))
# config.print_params(logger.info)


def main():
    logger.info("Logger is set - training start")

    # set gpu device id
    torch.cuda.set_device(config.gpu)
    
    # set seed
    np.random.seed(config.seed)
    torch.manual_seed(config.seed)
    torch.cuda.manual_seed_all(config.seed)

    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.enabled = True

    
    # get data with meta info
    input_size, input_channels, n_classes, train_data = utils.get_data(
        config.train_data, config.train_label, config.data_path, logger, cutout_length=0, validation=False)

    net_crit = nn.BCEWithLogitsLoss().to(device)
    model = SearchCNN(input_channels, config.init_channels, n_classes, config.layers, net_crit)

    try:
        logger.info("all gpus: {}".format(torch.cuda.device_count()))
        # model = nn.DataParallel(model, device_ids=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31] )
        # model = nn.DataParallel(model, device_ids=range(torch.cuda.device_count()))
        # model = nn.DataParallel(model)
        torch.distributed.init_process_group(backend="nccl")
        
        logger.info('1')
        model = model.to(device)
        logger.info('2')
        model = torch.nn.parallel.DistributedDataParallel(model, device_ids=range(torch.cuda.device_count()))
            
        
            

        # weights optimizer
        w_optim = torch.optim.SGD(model.module.weights(), config.w_lr, momentum=config.w_momentum,
                                  weight_decay=config.w_weight_decay)
        logger.info('3')

        # alphas optimizer
        alpha_optim = torch.optim.Adam(model.module.alphas(), config.alpha_lr, betas=(0.5, 0.999),
                                       weight_decay=config.alpha_weight_decay)
        logger.info('4')


        # split data to train/validation
        n_train = len(train_data)
        split = int(0.8 * n_train)
        indices = list(range(n_train))
        shuffle(indices)
        trainIndices = indices[:split]
        testIndices = indices[split:]
        with open(config.data_path+"trainTestIndices.pickle", "wb") as indicesFile:
            pickle.dump(trainIndices, indicesFile)
            pickle.dump(testIndices, indicesFile)

        with open(config.data_path+"trainTestIndices.pickle", "rb") as indicesFile:
            trainIndices = pickle.load(indicesFile)
        n_train = len(trainIndices)
        split = n_train // 2

        train_sampler = torch.utils.data.sampler.SubsetRandomSampler(trainIndices[:split])
        valid_sampler = torch.utils.data.sampler.SubsetRandomSampler(trainIndices[split:])
        train_loader = torch.utils.data.DataLoader(train_data,
                                                   batch_size=config.batch_size,
                                                   sampler=train_sampler,
                                                   num_workers=config.workers,
                                                   pin_memory=True)
        valid_loader = torch.utils.data.DataLoader(train_data,
                                                   batch_size=config.batch_size,
                                                   sampler=valid_sampler,
                                                   num_workers=config.workers,
                                                   pin_memory=True)
        logger.info('5')

        lr_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            w_optim, config.epochs, eta_min=config.w_lr_min)
        logger.info('6')

        architect = Architect(model, config.w_momentum, config.w_weight_decay)
        logger.info('7')

        # training loop
        best_genotype = None
        best_top1 = 0.
        for epoch in range(config.epochs):
            lr_scheduler.step()
            lr = lr_scheduler.get_lr()[0]

            model.module.print_alphas()
            logger.info('8')

            # training
            train(train_loader, valid_loader, model, architect, w_optim, alpha_optim, lr, epoch)
            logger.info('9')

            # validation
            cur_step = (epoch+1) * len(train_loader)
            top1 = validate(valid_loader, model, epoch, cur_step)
            logger.info('10')

            # log
            # genotype
            genotype = model.module.genotype()
            logger.info("genotype = {}".format(genotype))

            # genotype as a image
            plot_path = os.path.join(config.plot_path, "EP{:02d}".format(epoch+1))
            caption = "Epoch {}".format(epoch+1)
            plot(genotype.normal, plot_path + "-normal", caption)
            plot(genotype.reduce, plot_path + "-reduce", caption)

            # save
            if best_top1 < top1:
                best_top1 = top1
                best_genotype = genotype
                is_best = True
            else:
                is_best = False
            utils.save_checkpoint(model, config.path, is_best)
            print("")

        logger.info("Final best Prec@1 = {:.4%}".format(best_top1))
        logger.info("Best Genotype = {}".format(best_genotype))
    except Exception as e:
        logger.info("error: {}".format(e))



def train(train_loader, valid_loader, model, architect, w_optim, alpha_optim, lr, epoch):
    stats = utils.SumMeter()
    logger.info('8.1')
    # top1 = utils.AverageMeter()
    # top5 = utils.AverageMeter()
    losses = utils.AverageMeter()
    logger.info('8.2')

    cur_step = epoch*len(train_loader)
    # writer.add_scalar('train/lr', lr, cur_step)
    logger.info('8.3')

    model.train()
    logger.info('8.4')

    train_iter = iter(train_loader)
    valid_iter = iter(valid_loader)
    logger.info('8.5')
    
    for step, ((trn_X, trn_y), (val_X, val_y)) in enumerate(zip(train_loader, valid_loader)):
        if not torch.cuda.is_available():
            logging.info('no gpu device available')
            sys.exit(1)

        logger.info('8.6')

        trn_X, trn_y = trn_X.to(device, non_blocking=True), trn_y.to(device, non_blocking=True)
        logger.info('8.7')
        val_X, val_y = val_X.to(device, non_blocking=True), val_y.to(device, non_blocking=True)
        logger.info('8.8')
        N = trn_X.size(0)
        logger.info('8.9')
                        
        # phase 2. architect step (alpha)
        alpha_optim.zero_grad()
        logger.info('8.10')
        architect.unrolled_backward(trn_X, trn_y, val_X, val_y, lr, w_optim, logger)
        logger.info('8.11')
        alpha_optim.step()
        logger.info('8.12')

        # phase 1. child network step (w)
        w_optim.zero_grad()
        logger.info('8.13')
        logits = model(trn_X)
        logger.info('8.14')
        loss = model.module.criterion(logits, trn_y)
        logger.info('8.15')
        loss.backward()
        logger.info('8.16')
        # gradient clipping
        nn.utils.clip_grad_norm_(model.module.weights(), config.w_grad_clip)
        logger.info('8.17')
        w_optim.step()
        logger.info('8.18')

        truePos, trueNeg, falsePos, falseNeg = utils.accuracy(logits, trn_y)
        logger.info('8.19')
        losses.update(loss.item(), N)
        logger.info('8.20')
        stats.update(truePos, trueNeg, falsePos, falseNeg)
        logger.info('8.21')

        if step % config.print_freq == 0 or step == len(train_loader)-1:
            logger.info(
                "Train: [{:2d}/{}] Step {:03d}/{:03d} Loss {losses.avg:.3f} "
                "Prec@(1,5) ({top1:.1%}, {top5:.3})".format(
                    epoch+1, config.epochs, step, len(train_loader)-1, losses=losses,
                    top1=stats.accuracy(), top5=stats.MCC()))

        # writer.add_scalar('train/loss', loss.item(), cur_step)
        # writer.add_scalar('train/top1', prec1, cur_step)
        # writer.add_scalar('train/top5', prec5, cur_step)
        cur_step += 1

    logger.info("Train: [{:2d}/{}] Final Prec@1 {:.4}".format(epoch+1, config.epochs, stats.MCC()))


def validate(valid_loader, model, epoch, cur_step):
    stats = utils.SumMeter()
    # top1 = utils.AverageMeter()
    # top5 = utils.AverageMeter()
    losses = utils.AverageMeter()

    model.eval()

    with torch.no_grad():
        for step, (X, y) in enumerate(valid_loader):
            X, y = X.to(device, non_blocking=True), y.to(device, non_blocking=True)
            N = X.size(0)

            logits = model(X)
            loss = model.module.criterion(logits, y)

            truePos, trueNeg, falsePos, falseNeg = utils.accuracy(logits, y)
            losses.update(loss.item(), N)
            stats.update(truePos, trueNeg, falsePos, falseNeg)

            if step % config.print_freq == 0 or step == len(valid_loader)-1:
                logger.info(
                    "Valid: [{:2d}/{}] Step {:03d}/{:03d} Loss {losses.avg:.3f} "
                    "Prec@(1,5) ({top1:.1%}, {top5:.3})".format(
                        epoch+1, config.epochs, step, len(valid_loader)-1, losses=losses,
                        top1=stats.accuracy(), top5=stats.MCC()))

    # writer.add_scalar('val/loss', losses.avg, cur_step)
    # writer.add_scalar('val/top1', top1.avg, cur_step)
    # writer.add_scalar('val/top5', top5.avg, cur_step)

    logger.info("Valid: [{:2d}/{}] Final Prec@1 {:.4}".format(epoch+1, config.epochs, stats.MCC()))

    return stats.accuracy()


if __name__ == "__main__":
    main()
