#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=Bowtie2-Mapping
#SBATCH --array=0

source ~/.bashrc
source activate bowtie

# Define paths
job=/athena/masonlab/scratch/users/chb4004/PhD/bgc_datasets/monumentome/human_read_removal/job
names=($(cat $job))

BOWTIE2_INDEX=/athena/masonlab/scratch/users/chb4004/PhD/database/bowtie2_human/GRCh38_noalt_as/GRCh38_noalt_as
INPUT_DIR=/athena/masonlab/scratch/users/chb4004/PhD/bgc_datasets/monumentome/trim-gallore
OUTPUT_DIR=/athena/masonlab/scratch/users/chb4004/PhD/bgc_datasets/monumentome/human_read_removal

# Get sample name based on SLURMi task ID
SAMPLE_NAME=${names[${SLURM_ARRAY_TASK_ID}]}
R1=$INPUT_DIR/${SAMPLE_NAME}_L004_R1_001_val_1.fq.gz            #_1_val_1.fq.gz
R2=$INPUT_DIR/${SAMPLE_NAME}_L004_R2_001_val_2.fq.gz            #_2_val_2.fq.gz

ls -ld $R1 $R2

cd $OUTPUT_DIR

# Step 1: Bowtie2 Mapping (keeping both aligned and unaligned reads)
bowtie2 -p 8 -x $BOWTIE2_INDEX \
    -1 $R1 \
    -2 $R2 \
    -S $OUTPUT_DIR/${SAMPLE_NAME}_mapped_and_unmapped.sam

# Step 2: Convert SAM to BAM
samtools view -bS $OUTPUT_DIR/${SAMPLE_NAME}_mapped_and_unmapped.sam > $OUTPUT_DIR/${SAMPLE_NAME}_mapped_and_unmapped.bam
rm $OUTPUT_DIR/${SAMPLE_NAME}_mapped_and_unmapped.sam

# Step 3: Extract unmapped reads (both R1 and R2 unmapped)
samtools view -b -f 12 -F 256 $OUTPUT_DIR/${SAMPLE_NAME}_mapped_and_unmapped.bam > $OUTPUT_DIR/${SAMPLE_NAME}_bothReadsUnmapped.bam

# Step 4: Sort BAM file by read name
samtools sort -n -m 5G -@ 2 $OUTPUT_DIR/${SAMPLE_NAME}_bothReadsUnmapped.bam -o $OUTPUT_DIR/${SAMPLE_NAME}_bothReadsUnmapped_sorted.bam

# Step 5: Convert BAM back to paired-end FASTQ files
samtools fastq -@ 8 $OUTPUT_DIR/${SAMPLE_NAME}_bothReadsUnmapped_sorted.bam \
    -1 $OUTPUT_DIR/${SAMPLE_NAME}_host_removed_R1.fastq.gz \
    -2 $OUTPUT_DIR/${SAMPLE_NAME}_host_removed_R2.fastq.gz \
    -0 /dev/null -s /dev/null -n

echo Completed processing for $SAMPLE_NAME
