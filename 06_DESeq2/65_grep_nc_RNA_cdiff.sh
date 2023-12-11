#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=65_grep_nc_RNA_cdiff

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=fast
#SBATCH --array=1-3            # Array range times_series
##SBATCH --array=1-4             # Array range all except times_series
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
#rm *cdiff*/tables/*ncRNA*
######all except times_series######
#FILE_IN=$(ls sick*cdiff*/tables/*_reduced_final* early_late_cdiff/tables/*_reduced_final* IV_MI_cdiff/tables/*_reduced_final* | sed -n ${SLURM_ARRAY_TASK_ID}p)
#FILE_OUT=${FILE_IN/.complete_annot_reduced_final.tsv/_ncRNA.tsv}
######times_series######
FILE_IN=$(ls time_series_cdiff/tables/*.complete_annot_final.tsv | sed -n ${SLURM_ARRAY_TASK_ID}p)
FILE_OUT=${FILE_IN/.complete_annot_final.tsv/_ncRNA.tsv}
grep 'CD630_n\|CD630_s\|CD630_SQ\|CD630_Cdi' ${FILE_IN} > ${FILE_OUT}
