#!/bin/bash
#SBATCH --job-name=24_uninterleave_fastq
#SBATCH --partition=fast
#SBATCH --array=1-16
#SBATCH --cpus-per-task=8
#SBATCH --mem=2GB
#SBATCH --output=%x.%A_%a.stdout
#SBATCH --error=%x.%A_%a.stderr


# modules/dependencies loading

# Main script
READIN=$(ls *.cutadapt.trim.non_rrna.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

READOUT1=${READIN/.cutadapt.trim.non_rrna.fq.gz/.cleaned_R1.fq}
READOUT2=${READOUT1/_R1/_R2}

zcat -c $READIN | paste - - - - - - - - | tee >(cut -f1-4 | tr '\t' '\n' > ${READOUT1}) | cut -f5-8 | tr '\t' '\n' > ${READOUT2}

gzip ${READOUT1}
gzip ${READOUT2}

