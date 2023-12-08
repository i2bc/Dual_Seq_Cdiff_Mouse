#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=60_genesToExclude

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

# Objective: suppress tRNA & rRNA counts from counts files
# mouse : 486 genes, 64 rRNA + 422 tRNA 
# Clostri : 119 genes, 32 rRNA + 87 tRNA 

# based on "rRNA" and "tRNA" feature of gff files + Parent field (for gene ID)
# create a file of one line (where each gene name is coma separated)
srun awk -F "\t" '{if(($3=="tRNA")||($3=="rRNA")){print $0}}' ../04_featureCounts_genome/mergeAnnot.gff | grep -o "Parent=[^;]*;" | sed 's/Parent=//;s/;//' | awk 'BEGIN{RS="\n";ORS=","}{print}' | sed 's/,$/\n/' > genesToExclude.txt


