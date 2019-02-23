#!/bin/bash
#PBS -l wd
#PBS -o deepBrain_train.out
#PBS -e deepBrain_train.err
#PBS -q gpu
#PBS -l mem=256GB
#PBS -l ngpus=6
#PBS -l ncpus=18
#PBS -P yr31
#PBS -l walltime=48:00:00

module load pytorch

python3 /short/yr31/aa7970/azData/DeepBrain/Scripts/search.py --name deepbrain --datapath Data/ --train_data H3K27ac_rnaSeq.Pos.Neg.tfSpecificInput.memmap --train_label H3K27ac_rnaSeq.Pos.Neg.tfSpecificLabel.memmap --epochs 5 --batch_size 32