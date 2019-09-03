#!/bin/sh

#PBS -N Split_Numpy_katana
#PBS -o Split_Numpy.out
#PBS -e Split_Numpy.err
#PBS -l select=ncpus=8:mem=190G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae

module load python/3.7.3

cd $PBS_O_WORKDIR

python3 /srv/scratch/z3526914/DeepBrain/Scripts/split_and_Numpy_V2.py	\
    --datadir /srv/scratch/z3526914/DeepBrain/Data/  \
    --datafilename HumanFC_ENCODE_EpiMap_nonZero.bin.Seq_Labels  \
    --valid_chr_id 1
