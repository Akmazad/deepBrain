#!/bin/bash
#PBS -l wd
#PBS -o deepBrain_Split_and_Numpy.out
#PBS -e deepBrain_Split_and_Numpy.err
#PBS -q hugemem
#PBS -l mem=256GB
#PBS -l ncpus=14
#PBS -P yr31
#PBS -l walltime=1:00:00


#module load python3/3.7.2
#export PYTHONPATH=/home/561/aa7970/.local/lib/python3.7/site-packages:$PYTHONPATH

#module unload gcc
#module unload intel-fc intel-cc intel-mkl
#
#
module load gcc/6.2.0
module load python3/3.6.7-gcc620
export LD_PRELOAD=/apps/gcc/6.2.0/lib64/libstdc++.so.6
export PYTHONPATH=/short/yr31/aa7970/local/lib/python3.6/site-packages:$PYTHONPATH


python3 /short/yr31/aa7970/azData/DeepBrain/Scripts/split_and_Numpy.py  \
    --datadir /short/yr31/aa7970/azData/DeepBrain/Data/  \
    --datafilename H3K27ac_rnaSeq.Pos.tfSpecific  \
    --valid_chr_id 1