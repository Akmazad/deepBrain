#!/bin/sh
#PBS -l wd
#PBS -o humanFC.out
#PBS -e humanFC.err
#PBS -q hugemem
#PBS -l mem=1TB
#PBS -l ncpus=7
#PBS -P yr31
#PBS -l walltime=06:00:00

module load R/3.5.1

Rscript /short/yr31/aa7970/azData/DeepBrain/Scripts/HumanFC_post_processing.R