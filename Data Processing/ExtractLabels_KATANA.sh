#!/bin/sh

#PBS -N ExtractLabels_katana
#PBS -o ExtractLabels.out
#PBS -e ExtractLabels.err
#PBS -l select=ncpus=8:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae


module load R/3.5.3

Rscript /srv/scratch/z3526914/DeepBrain/Scripts/ExtractLabels_KATANA.R \
	/srv/scratch/z3526914/DeepBrain/Data/ \
	HumanFC_tf_specific_labels.bed \
	EpiMap_tf_specific_labels.bed \
	CAGE_tf_specific_labels.bed \
	ENCODE_TFs_tf_specific_labels.bed \
	HumanFC_CAGE_ENCODE_EpiMap_tf_specific.bin.Seq.bed \
	HumanFC_CAGE_ENCODE_EpiMap_tf_specific.bin.Labels.bed \
	HumanFC_CAGE_ENCODE_EpiMap_tf_specific.bin.Seq_Labels.bed
