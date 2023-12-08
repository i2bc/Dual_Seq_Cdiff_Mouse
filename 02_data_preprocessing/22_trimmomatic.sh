#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=22_trimmomatic

### Requirements
#SBATCH --partition=fast
#SBATCH --array=1-16
#SBATCH --cpus-per-task=8
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
module load trimmomatic


# Main scrip
R1IN=$(ls *_R1.cutadapt.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
R2IN=${R1IN/_R1/_R2}
BASEOUT=${R1IN/_R1.cutadapt.fq.gz/.cutadapt.trim.fq.gz}
STATOUT=${BASEOUT/.fq.gz/.stats}

echo $R1IN
echo $R2IN
echo $BASEOUT
echo $STATOUT

srun trimmomatic PE -threads ${SLURM_CPUS_PER_TASK} -summary $STATOUT -validatePairs $R1IN $R2IN -baseout $BASEOUT LEADING:20 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:10


echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
