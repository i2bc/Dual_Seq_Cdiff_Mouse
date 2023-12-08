#!/bin/bash
#SBATCH --job-name=33_bowtie2_db_build
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --output=%x.%j.stdout    # Standard output for single jobs
#SBATCH --error=%x.%j.stderr     # Standard error for single jobs

# modules/dependencies loading
#  bowtie2/2.4.1
module load bowtie2

# Main script

srun bowtie2-build --threads ${SLURM_CPUS_PER_TASK} Mmusculus_GRCm39_Cdifficile_v2_merged.fa Mmusculus_GRCm39_Cdifficile_v2_merged.fa.bowtie2_index
