#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=64_modification_table_text

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=fast
#SBATCH --array=1-50                # Array range
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB

### Email
##SBATCH --mail-user=email@address
##SBATCH --mail-type=ALL

### Output
##SBATCH --output=%x.%j.stdout    # Standard output for single jobs
##SBATCH --error=%x.%j.stderr     # Standard error for single jobs
#SBATCH --output=%x.%A_%a.stdout    # Standard output for array jobs
#SBATCH --error=%x.%A_%a.stderr     # Standard error for array jobs


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

# usage example:

# Main script
#rm */tables/*final*
#rm *tables/*tmp*
FILE_IN=$(ls */tables/*annot*tsv | sed -n ${SLURM_ARRAY_TASK_ID}p)
FILE_TMP=${FILE_IN/.tsv/.tmp}
srun sed "s/\./\,/g;s/%27/'/g;s/%28/(/g;s/%29/)/g;s/%2B/+/g;s/%2C/,/g;s/%2D/-/g;s/%3B/;/g;s/%5B/[/g;s/%5D/]/g" ${FILE_IN} > ${FILE_TMP}
#rename .tmp _final.tsv */tables/*tmp
