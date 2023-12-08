#!/bin/bash
#SBATCH --job-name=31_hisat2_db_build
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --output=%x.%j.stdout    # Standard output for single jobs
#SBATCH --error=%x.%j.stderr     # Standard error for single jobs

# modules/dependencies loading
module load hisat2

# Main script
srun hisat2-build -p ${SLURM_CPUS_PER_TASK} Mmusculus_GRCm39_Cdifficile_v2_merged.fa Mmusculus_GRCm39_Cdifficile_v2_merged.hisat2_index
