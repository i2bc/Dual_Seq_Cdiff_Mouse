#!/bin/bash
#SBATCH --job-name=32_hisat2_mapping
#SBATCH --partition=fast
#SBATCH --array=1-16
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --output=%x.%A_%a.stdout
#SBATCH --error=%x.%A_%a.stderr

# modules/dependencies loading
module load hisat2 samtools

# Main script
R1_in=$(ls *.cleaned_R1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
R2_in=${R1_in/_R1/_R2}

basename_out=${R1_in/.cleaned_R1.fq.gz/.cleaned.hisat2}
basename_unmap=${R1_in/.cleaned_R1.fq.gz/.unmaped_R%.hisat2}

echo $R1_in
echo $R2_in
echo $basename_out
echo $basename_unmap

srun hisat2 -p ${SLURM_CPUS_PER_TASK} -x Mmusculus_GRCm39_Cdifficile_v2_merged.hisat2_index --un-conc-gz ${basename_unmap}.gz -1 $R1_in -2 $R2_in | tee >(samtools flagstat - > ${basename_out}.flagstat) | samtools sort -O BAM > ${basename_out}.bam

srun samtools index ${basename_out}.bam
