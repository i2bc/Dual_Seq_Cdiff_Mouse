#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=50_S9_samtools_merge

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=fast
##SBATCH --array=1-2                # Array range
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=2GB

### Email
##SBATCH --mail-user=email@address
##SBATCH --mail-type=ALL

### Output
#SBATCH --output=%x.%j.stdout    # Standard output for single jobs
#SBATCH --error=%x.%j.stderr     # Standard error for single jobs
##SBATCH --output=%x.%A_%a.stdout    # Standard output for array jobs
##SBATCH --error=%x.%A_%a.stderr     # Standard error for array jobs

################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'

# modules/dependencies loading
#module load subread/2.0.1
module load samtools

# Main script

srun samtools merge S9_S46.cleaned.hisat2_sorted.bam S9_S4.cleaned.hisat2_sorted.bam S9_S6.cleaned.hisat2_sorted.bam
srun samtools merge S9_S46.hisat2_unmaped.bowtie2_sorted.bam S9_S4.hisat2_unmaped.bowtie2_sorted.bam S9_S6.hisat2_unmaped.bowtie2_sorted.bam
