#!/bin/bash
#PBS -l wd
#PBS -o deepBrain_train.out
#PBS -e deepBrain_train.err
#PBS -q gpu
#PBS -l mem=100GB
#PBS -l ngpus=8
#PBS -l ncpus=48
#PBS -P yr31
#PBS -l walltime=48:00:00

module unload gcc
module unload python3/3.6.2
module unload intel-fc intel-cc intel-mkl

module load gcc/6.2.0
module load python3/3.6.7-gcc620
export PYTHONPATH=/short/yr31/aa7970/local/lib/python3.6/site-packages:$PYTHONPATH
export LD_PRELOAD=/apps/gcc/6.2.0/lib64/libstdc++.so.6


python3 /short/yr31/aa7970/azData/DeepBrain/Scripts/search.py --name deepbrain --datapath Data/ --train_data H3K27ac_rnaSeq.Pos.Neg.tfSpecificInput.memmap --train_label H3K27ac_rnaSeq.Pos.Neg.tfSpecificLabel.memmap --epochs 5 --batch_size 64