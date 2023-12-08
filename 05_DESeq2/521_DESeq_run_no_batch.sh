#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=621_DESeq_run_no_batch

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=fast
##SBATCH --array=1-10                # Array range
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB

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
module load r/3.6.3
module load pandoc/2.7.2

# Main script
# $1 : design file (to create result directory)
# $2 : value of --condRef
# $3 : no batch effect allowed (NULL)
# usage : sbatch script.sh $1 $2 $3
# exemple : sbatch 621_DESeq_run.sh allSamples_mouse_counts.txt MI32h NULL

analyse=`basename $1 .txt`
rm -Rf ${analyse}
mkdir -p ${analyse}
cd ${analyse}

Rscript ../sartools_script_DESeq2_CL.r --projectName=${analyse} --targetFile=../${analyse}.txt --rawDir=../fc_data --varInt=group --condRef=$2
#mv tables figures ${analyse}.RData ${analyse}_report.html ${analyse}/.

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
