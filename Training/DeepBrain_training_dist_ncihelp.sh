#!/bin/bash -l
#PBS -l wd
#PBS -o deepBrain_train_dist_hcihelp.out
#PBS -e deepBrain_train_dist_hcihelp.err
#PBS -l ncpus=48
#PBS -l ngpus=16
#PBS -P yr31
#PBS -l walltime=01:30:00
#PBS -l mem=50GB
#PBS -q gpu

# Get number of CPUs per GPU.
# This is used as the "step" when indexing the nodes with pbsdsh so that we get the right number per node.

cpus_per_gpu=$((PBS_NCPUS / PBS_NGPUS))
master_addr=$(grep $(hostname --short) /etc/hosts | awk '{print $1}')
port=29500
nNode=2

# Launch the processes on all nodes and wait for it to complete.
for rank in $(seq 0 $((nNode - 1))); do
  pbsdsh -n $((rank * 24)) -- /short/yr31/aa7970/azData/DeepBrain/Scripts/launch.sh ${master_addr} ${port} ${nNode} ${rank} &
done
wait