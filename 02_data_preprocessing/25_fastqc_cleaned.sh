#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=25_fastqc_cleaned

### Requirements
#SBATCH --partition=fast
#SBATCH --array=1-32
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB

### Output
#SBATCH --output=%x.%A_%a.stdout
#SBATCH --error=%x.%A_%a.stderr

################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'

# modules loading
module load fastqc


# Main script
READIN=$(ls *cleaned_R1.fq.gz *cleaned_R2.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
srun fastqc $READIN


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
