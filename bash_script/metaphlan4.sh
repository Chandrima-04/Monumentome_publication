#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --job-name=Metaphlan_4
#SBATCH --array=0-231

source ~/.bashrc
source activate metaphlan4

job=/athena/masonlab/scratch/users/chb4004/PhD/monumentome/job
names=($(cat $job))

input1=/athena/masonlab/scratch/users/chb4004/PhD/bgc_datasets/monumentome/human_read_removal/${names[${SLURM_ARRAY_TASK_ID}]}_host_removed_R1.fastq.gz
input2=/athena/masonlab/scratch/users/chb4004/PhD/bgc_datasets/monumentome/human_read_removal/${names[${SLURM_ARRAY_TASK_ID}]}_host_removed_R2.fastq.gz

ls -ld $input1 $input2

cd /athena/masonlab/scratch/users/chb4004/PhD/bgc_datasets/monumentome/metaphlan4_17082026_alltaxa
work_dir=${names[${SLURM_ARRAY_TASK_ID}]}
rm -rf $work_dir
mkdir -p $work_dir
cd $work_dir

merged_fastq_gz=merged-${SLURM_ARRAY_TASK_ID}.fastq.gz

zcat $input1 $input2 | gzip -c > $merged_fastq_gz

metaphlan --input_type fastq --bowtie2db /athena/masonlab/scratch/users/chb4004/databases/metaphlan_vJan25/    --nproc 10  $merged_fastq_gz ${names[${SLURM_ARRAY_TASK_ID}]}.tsv
