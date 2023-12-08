#!/bin/bash
#SBATCH --job-name=34_bowtie2_mapping
#SBATCH --partition=fast
#SBATCH --array=1-16
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --output=%x.%A_%a.stdout
#SBATCH --error=%x.%A_%a.stderr

# modules/dependencies loading
# bowtie2/2.4.1
# samtools/1.10
module load bowtie2 samtools

# Main script
R1_in=$(ls *.unmaped_R1.hisat2.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
R2_in=${R1_in/_R1/_R2}

basename_out=${R1_in/.unmaped_R1.hisat2.gz/.hisat2_unmaped.bowtie2}
basename_unmap=${R1_in/.unmaped_R1.hisat2.gz/.unmaped_R%.bowtie2}

srun bowtie2 -p ${SLURM_CPUS_PER_TASK} --omit-sec-seq --sam-no-qname-trunc --un-conc-gz ${basename_unmap}.gz -x Mmusculus_GRCm39_Cdifficile_v2_merged.fa.bowtie2_index -1 $R1_in -2 $R2_in | tee >(samtools flagstat - > ${basename_out}.flagstat) | samtools sort -O BAM > ${basename_out}.bam

srun samtools index ${basename_out}.bam
