# deep Brain

## Step-1: Preparing RNA-seq data from Brain tissue for Deep learning applications

### Code File: 
1) "windowMat.r" contains all the functions required to RNA-seq preparation (note: this code should be scaled up for running on raijin, including reading files from remote host()
2) "sample_code1.r" contains very basic codes demonstrating "derfinder" mainly

## Step-2: 
PBS script for running on HPC
#!/bin/bash
#PBS -P yr31
#PBS -q gpu
#PBS -l ngpus=2
#PBS -l ncpus=6
#PBS -l walltime=0:45:00,mem=8GB
#PBS -l wd

R rScript.r
