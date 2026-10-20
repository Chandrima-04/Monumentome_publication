#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=Trim-Galore
#SBATCH --array=0-231
#SBATCH --mail-user=cb4603@nyu.edu
#SBATCH --mail-type=FAIL

source ~/.bashrc
source activate trim_gallore

job=/athena/masonlab/scratch/users/chb4004/PhD/monumentome/job
names=($(cat $job))

input1=/athena/masonlab/scratch/users/chb4004/PhD/monumentome/data/${names[${SLURM_ARRAY_TASK_ID}]}_L004_R1_001.fastq.gz
input2=/athena/masonlab/scratch/users/chb4004/PhD/monumentome/data/${names[${SLURM_ARRAY_TASK_ID}]}_L004_R2_001.fastq.gz

ls -ld $input1 $input2

cd /athena/masonlab/scratch/users/chb4004/PhD/bgc_datasets/monumentome/trim-gallore

trim_galore --paired $input1 $input2
