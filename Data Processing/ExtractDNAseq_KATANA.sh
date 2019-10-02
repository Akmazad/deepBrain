#!/bin/sh

#PBS -N ExtractDNA_katana
#PBS -o ExtractDNA.out
#PBS -e ExtractDNA.err
#PBS -l select=ncpus=8:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae


module load R/3.5.3

Rscript /srv/scratch/z3526914/DeepBrain/Scripts/ExtractDNAseq_KATANA.R \
	/srv/scratch/z3526914/DeepBrain/Data \
	400 \
	HumanFC_ENCODE_EpiMap_nonZero.binInfo.Union.bed \
	HumanFC_ENCODE_EpiMap_nonZero.bin.Seq.bed
