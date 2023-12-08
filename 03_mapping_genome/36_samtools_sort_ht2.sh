#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=36_samtools_sort_ht2

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=fast
#SBATCH --array=1-16                # Array range
##SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50GB

### Email
##SBATCH --mail-user=email@address
##SBATCH --mail-type=ALL

### Output
##SBATCH --output=%x.%j.stdout    # Standard output for single jobs
##SBATCH --error=%x.%j.stderr     # Standard error for single jobs
#SBATCH --output=%x.%A_%a.stdout    # Standard output for array jobs
#SBATCH --error=%x.%A_%a.stderr     # Standard error for array jobs

### Workdir
###SBATCH --workdir=/shared/projects/intact-lt/3_count_features/

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
module load samtools

# Main script

in_file=$(ls 32_hisat2_mapping_bam/*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_file=${in_file/.bam/_sorted.bam}
# sorted by read name (-n)
srun samtools sort -n -@ ${SLURM_CPUS_PER_TASK} -o ${out_file} ${in_file}
# Déplacement manuel des _sorted.bam dans le dossier 36_samtools_sorted_ht2_bam/

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)

