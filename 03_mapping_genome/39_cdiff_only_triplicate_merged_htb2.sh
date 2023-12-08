#!/bin/bash
#SBATCH --job-name=cdiff_only_triplicate_merged_htb2
#SBATCH --partition=fast
#SBATCH --array=1-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --output=%x.%A_%a.stdout
#SBATCH --error=%x.%A_%a.stderr

# modules/dependencies loading
# samtools/1.10
module load samtools

# Main script

FILE_IN=$(ls 39_cdiff_only_triplicate_merged_htb2/*.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
FILE_OUT_TMP=${FILE_IN/.txt/.tmp.bam}
FILE_OUT=${FILE_IN/.txt/.bam}
srun samtools merge -h 39_cdiff_only_triplicate_merged_htb2/header.sam -R NC_009089.1 -b ${FILE_IN} ${FILE_OUT_TMP}
srun samtools view -F 0x100 -b -o ${FILE_OUT} ${FILE_OUT_TMP}
srun samtools index ${FILE_OUT}
rm ${FILE_OUT_TMP}
