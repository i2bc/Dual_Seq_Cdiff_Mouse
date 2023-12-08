#!/bin/bash
#SBATCH --job-name=IV_sampled_bam
#SBATCH --partition=fast
#SBATCH --array=1-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --output=%x.%A_%a.stdout
#SBATCH --error=%x.%A_%a.stderr

# modules/dependencies loading
# samtools/1.10
module load samtools

# Main script

READIN_ht2=$(ls 32_hisat2_mapping_bam/IV*.cleaned.hisat2.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
READIN_bt2=$(ls 34_bowtie2_mapping_bam/IV*.hisat2_unmaped.bowtie2.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
SAMPLE_NAME_ht2=$(echo $READIN_ht2|sed 's/.cleaned.hisat2.bam/_sampled.cleaned.hisat2.bam/;s/32_hisat2_mapping_bam/38_IV_sampled_bam/')
SAMPLE_NAME_bt2=$(echo $READIN_bt2|sed 's/.hisat2_unmaped.bowtie2.bam/_sampled.hisat2_unmaped.bowtie2.bam/;s/34_bowtie2_mapping_bam/38_IV_sampled_bam/')
OUT_FILE_ht2=$(echo ${SAMPLE_NAME_ht2}|sed 's/\.bam/_sorted.bam/')
OUT_FILE_bt2=$(echo ${SAMPLE_NAME_bt2}|sed 's/\.bam/_sorted.bam/')
OUT_FILE_ht2_pos=$(echo ${SAMPLE_NAME_ht2}|sed 's/\.bam/_sorted_pos.bam/')
OUT_FILE_bt2_pos=$(echo ${SAMPLE_NAME_bt2}|sed 's/\.bam/_sorted_pos.bam/')

# Sampling IV
#echo ${READIN_ht2} ${SAMPLE_NAME_bt2} ${OUT_FILE_ht2} 
#echo ${READIN_bt2} ${SAMPLE_NAME_bt2} ${OUT_FILE_bt2}
#srun samtools view -F 0x100 -s 1.005 -b ${READIN_ht2} > ${SAMPLE_NAME_ht2}
#srun samtools view -F 0x100 -s 1.005 -b ${READIN_bt2} > ${SAMPLE_NAME_bt2}

# Sorting IV_sampled by name (-n)
#srun samtools sort -n -@ ${SLURM_CPUS_PER_TASK} -o ${OUT_FILE_ht2} ${SAMPLE_NAME_ht2}
#srun samtools sort -n -@ ${SLURM_CPUS_PER_TASK} -o ${OUT_FILE_bt2} ${SAMPLE_NAME_bt2}

# Sorting IV_sampled by position
srun samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${OUT_FILE_ht2_pos} ${SAMPLE_NAME_ht2}
srun samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${OUT_FILE_bt2_pos} ${SAMPLE_NAME_bt2}
srun samtools index ${OUT_FILE_ht2_pos} 
srun samtools index ${OUT_FILE_bt2_pos}
