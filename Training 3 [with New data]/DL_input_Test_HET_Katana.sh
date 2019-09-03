#!/bin/bash

#PBS -N HET_training_katana
#PBS -o deepBrain_katana.out
#PBS -e deepBrain_katana.err
#PBS -l select=ncpus=32:ngpus=4:mem=180G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae


module load pytorch


cd $PBS_O_WORKDIR

python3 /srv/scratch/z3526914/DeepBrain/Scripts/DL_input_Test_HET_Katana.py
