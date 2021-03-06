""" Config class for search/augment """
import argparse
import os
import genotypes as gt


class BaseConfig(argparse.Namespace):
    def print_params(self, prtf=print):
        prtf("")
        prtf("Parameters:")
        for attr, value in sorted(vars(self).items()):
            prtf("{}={}".format(attr.upper(), value))
        prtf("")

    def as_markdown(self):
        """ Return configs as markdown format """
        text = "|name|value|  \n|-|-|  \n"
        for attr, value in sorted(vars(self).items()):
            text += "|{}|{}|  \n".format(attr, value)

        return text


class SearchConfig(BaseConfig):
    def build_parser(self):
        parser = argparse.ArgumentParser("Search config")
        parser.add_argument('--name', required=True)
        parser.add_argument('--datapath', required=True)
        parser.add_argument('--train_data', required=True)
        parser.add_argument('--train_label', required=True)
        parser.add_argument('--w_lr', type=float, default=0.01, help='lr for weights')
        parser.add_argument('--w_lr_min', type=float, default=0.001, help='minimum lr for weights')
        parser.add_argument('--w_momentum', type=float, default=0.9)
        parser.add_argument('--w_weight_decay', type=float, default=3e-4)
        parser.add_argument('--w_grad_clip', type=float, default=5.,
                            help='gradient clipping for weights')
        parser.add_argument('--print_freq', type=int, default=50, help='print frequency')
        # parser.add_argument('--gpu', type=int, default=0, help='gpu device id')
        # parser.add_argument('--epochs', type=int, default=50, help='# of training epochs')
        parser.add_argument('--init_channels', type=int, default=16)
        parser.add_argument('--layers', type=int, default=8)
        parser.add_argument('--batch_size', type=int, default=256, metavar='N', help='mini-batch size (default: 256), this is the total batch size of all GPUs on the current node when using Data Parallel or Distributed Data Parallel')
        parser.add_argument('--seed', type=int, default=None, help='seed for initializing training')
        parser.add_argument('--workers', type=int, default=4, help='# of workers')
        parser.add_argument('--alpha_lr', type=float, default=3e-4, help='lr for alpha')
        parser.add_argument('--alpha_weight_decay', type=float, default=1e-4, help='weight decay for alpha')
        parser.add_argument('--world-size', default=-1, type=int, help='number of nodes for distributed training')
        parser.add_argument('--rank', default=-1, type=int, help='node rank for distributed training')
        parser.add_argument('--dist-url', default='env://', type=str, help='url used to set up distributed training')
        parser.add_argument('--dist-backend', default='nccl', type=str, help='distributed backend')
        parser.add_argument('--gpu', default=None, type=int, help='GPU id to use.')
        parser.add_argument('--multiprocessing-distributed', default=True, action='store_true', help='Use multi-processing distributed training to launch N processes per node, which has N GPUs. This is the fastest way to use PyTorch for either single node or multi node data parallel training')
        parser.add_argument('--pretrained', dest='pretrained', action='store_true', help='use pre-trained model')
        parser.add_argument('--epochs', default=90, type=int, metavar='N', help='number of total epochs to run')
        parser.add_argument('--start-epoch', default=0, type=int, metavar='N', help='manual epoch number (useful on restarts)')
        parser.add_argument('--resume', default='', type=str, metavar='PATH', help='path to latest checkpoint (default: none)')


        return parser

    def __init__(self):
        parser = self.build_parser()
        args = parser.parse_args()
        super().__init__(**vars(args))

        # first get the grand-parent directory of the current main-script (search.py) '
        self.data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # then join the datapath
        self.data_path = os.path.join(self.data_path, args.datapath)
        self.path = os.path.join(self.data_path, self.name, 'searchs')
        self.plot_path = os.path.join(self.path, 'plots')


class AugmentConfig(BaseConfig):
    def build_parser(self):
        parser = argparse.ArgumentParser("Augment config")
        parser.add_argument('--name', required=True)
        parser.add_argument('--datapath', required=True)
        parser.add_argument('--train_data', required=True)
        parser.add_argument('--train_label', required=True)
        parser.add_argument('--batch_size', type=int, default=96)
        parser.add_argument('--lr', type=float, default=0.025, help='lr for weights')
        parser.add_argument('--momentum', type=float, default=0.9)
        parser.add_argument('--weight_decay', type=float, default=3e-4)
        parser.add_argument('--grad_clip', type=float, default=5.,
                            help='gradient clipping for weights')
        parser.add_argument('--print_freq', type=int, default=200, help='print frequency')
        parser.add_argument('--gpu', type=int, default=None, help='gpu device id')
        parser.add_argument('--epochs', type=int, default=600, help='# of training epochs')
        parser.add_argument('--init_channels', type=int, default=36)
        parser.add_argument('--layers', type=int, default=20)
        parser.add_argument('--seed', type=int, default=None, help='random seed')
        parser.add_argument('--workers', type=int, default=4, help='# of workers')
        parser.add_argument('--aux_weight', type=float, default=0.4, help='auxiliary loss weight')
        parser.add_argument('--cutout_length', type=int, default=16, help='cutout length')
        parser.add_argument('--drop_path_prob', type=float, default=0.2, help='drop path prob')

        parser.add_argument('--genotype', required=True, help='Cell genotype')

        return parser

    def __init__(self):
        parser = self.build_parser()
        args = parser.parse_args()
        super().__init__(**vars(args))

        # first get the grand-parent directory of the current main-script (search.py) '
        self.data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # then join the datapath
        self.data_path = os.path.join(self.data_path, args.datapath)
        self.path = os.path.join(self.data_path, self.name, 'augments')
        self.genotype = gt.from_str(self.genotype)
