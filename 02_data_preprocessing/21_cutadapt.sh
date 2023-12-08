#!/bin/bash
#SBATCH --job-name=21_cutadapt
#SBATCH --partition=fast
#SBATCH --array=1-16
#SBATCH --cpus-per-task=8
#SBATCH --mem=12GB
#SBATCH --output=%x.%A_%a.stdout
#SBATCH --error=%x.%A_%a.stderr


module load cutadapt
R1_in=$(ls *_R1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
R2_in=${R1_in/_R1/_R2}

R1_out=${R1_in/_R1.fastq/_R1.cutadapt.fq}
R2_out=${R2_in/_R2.fastq/_R2.cutadapt.fq}

srun cutadapt -j ${SLURM_CPUS_PER_TASK} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $R1_out -p $R2_out $R1_in $R2_in
