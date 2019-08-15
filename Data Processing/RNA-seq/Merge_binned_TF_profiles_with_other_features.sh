#!/bin/sh
#PBS -l wd
#PBS -o mergeTFs.out
#PBS -e mergeTFs.err
#PBS -q hugemem
#PBS -l mem=1TB
#PBS -l ncpus=7
#PBS -P yr31
#PBS -l walltime=96:00:00

module load R/3.5.1

Rscript /short/yr31/aa7970/azData/DeepBrain/Data/Merge_binned_TF_profiles_with_other_features.R