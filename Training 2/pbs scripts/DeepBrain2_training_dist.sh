#!/bin/bash -l

#PBS -l wd
#PBS -o deepBrain2_train_dist_np.out
#PBS -e deepBrain2_train_dist_np.err
#PBS -l ncpus=24
#PBS -l ngpus=8
#PBS -P yr31
#PBS -l walltime=5:00:00
#PBS -l mem=256GB
#PBS -q gpu


module load python3/3.6.2
module load pytorch/1.1.0a0-py36
module load openmpi/3.1.3
module load cuda/9.0
module load cudnn/7.1.1-cuda9.0
module load gcc/4.9.0
module load cmake/3.8.2
module load magma/2.3.0
module load intel-mkl/17.0.1.132
module load hdf5/1.10.2
module load nccl/2.2.13-cuda9.0

# start
module unload gcc
module unload python3/3.6.2
module unload intel-fc intel-cc intel-mkl

module load gcc/6.2.0
module load python3/3.6.7-gcc620
export PYTHONPATH=/short/yr31/aa7970/local/lib/python3.6/site-packages:$PYTHONPATH
export LD_PRELOAD=/apps/gcc/6.2.0/lib64/libstdc++.so.6
# end

# Parse command line.
export MASTER_ADDR=$(grep $(hostname --short) /etc/hosts | awk '{print $1}')
export MASTER_PORT=29500
export WORLD_SIZE=1
export RANK=0

# Setup local details.
gpus_per_node=8
gpus_per_numa=4
gpu_id=$((RANK % gpus_per_node))
numa_node=$((gpu_id / gpus_per_numa))

# Setup other details.
export NCCL_DEBUG=INFO
export NCCL_MIN_NRINGS=5


# Launch program.
# exec numactl --membind=${numa_node} --cpunodebind=${numa_node} -- \
  python3 /short/yr31/aa7970/azData/DeepBrain/Scripts/deepbrain2_dist.py  \
    --name deepbrainStaticConvnet_chr1 \
    --DataDir Data/ \
    --TrainingDataFile H3K27ac_rnaSeq.Pos.tfSpecific_trainingData_chr1_value.npy \
    --TrainingLabelFile H3K27ac_rnaSeq.Pos.tfSpecific_trainingData_chr1_label_all.npy \
    --TestingDataFile H3K27ac_rnaSeq.Pos.tfSpecific_validationData_chr1_value.npy \
    --TestingLabelFile H3K27ac_rnaSeq.Pos.tfSpecific_validationData_chr1_label_all.npy \
    --nEpochs 1 \
    --BATCH_SIZE 128 \
    --NUM_OUTPUTS 131 \
    --multiprocessing-distributed \
