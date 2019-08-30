#!/bin/sh

#PBS -N ExtractLabels_katana
#PBS -o ExtractLabels.out
#PBS -e ExtractLabels.err
#PBS -l select=ncpus=8:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae


module load R/3.5.3

Rscript /srv/scratch/z3526914/DeepBrain/Scripts/ExtractLabels_KATANA.R