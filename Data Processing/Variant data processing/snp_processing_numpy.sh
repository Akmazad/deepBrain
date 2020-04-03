#!/bin/sh

#PBS -N snp_processing_numpy_katana
#PBS -o snp_processing_numpy.out
#PBS -e snp_processing_numpy.err
#PBS -l select=ncpus=8:mem=250G
#PBS -l walltime=6:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae

module load python/3.7.3

cd $PBS_O_WORKDIR

python3 /srv/scratch/z3526914/DeepBrain/Scripts/snp_processing_numpy.py	\
    --datadir /srv/scratch/z3526914/DeepBrain/Data/  \
    --datafilename deepsea_supple_tabale5_fastaseq
