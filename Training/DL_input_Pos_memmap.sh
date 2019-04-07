#!/bin/bash
#PBS -l wd
#PBS -o deepBrain_Pos_memMap.out
#PBS -e deepBrain_Pos_memMap.err
#PBS -q hugemem
#PBS -l mem=256GB
#PBS -l ncpus=7
#PBS -P yr31
#PBS -l walltime=4:00:00

module load pytorch

python3 /short/yr31/aa7970/azData/DeepBrain/Scripts/csv2memmap.py /short/yr31/aa7970/azData/DeepBrain/Data/ H3K27ac_rnaSeq.Pos.tfSpecific