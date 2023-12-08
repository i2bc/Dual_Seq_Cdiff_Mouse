#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=51_bht2_featureCounts

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=fast
##SBATCH --array=1-20                # Array range
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80GB

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
module load subread/2.0.1

# merge annotations files for the 2 genomes 
if [ ! -f mergeAnnot.gff ] 
then
   cat < /shared/projects/clos_ribo_reg/reference_dbs/Mmusculus/GRCm39/GCF_000001635.27_GRCm39_genomic.gff > mergeAnnot.gff 
   cat < /shared/projects/clos_ribo_reg/reference_dbs/Cdifficile_MaGe/CD630_MaGe_ID.gff >> mergeAnnot.gff
   # cat < plasmid >> mergeAnnot.gff
fi
# Main script
# counts from hisat2 mapping (dual species)
srun featureCounts -T ${SLURM_CPUS_PER_TASK} -p -s 2 -O -M --fraction --primary -t gene -g ID -a mergeAnnot.gff -o all_ht2_featureCounts.txt *.cleaned.hisat2_sorted.bam
# srun head -n 2 all_ht2_featureCounts.txt | tail -n 1 | cut -d $'\t' -f 1,7-27 | sed 's/.cleaned.hisat2_sorted.bam//g' > all_ht2_mouse_featureCounts_matrix.tab
# srun head -n 2 all_ht2_featureCounts.txt | tail -n 1 | cut -d $'\t' -f 1,7-27 | sed 's/.cleaned.hisat2_sorted.bam//g' > all_ht2_cdiff_featureCounts_matrix.tab
# srun awk -F "\t" '{if(NR>2){if(($2!="NC_008226.2")&&($2!="NC_009089.1")){print $1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27>>"all_ht2_mouse_featureCounts_matrix.tab"}else{print $1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27>>"all_ht2_cdiff_featureCounts_matrix.tab"}}}' all_ht2_featureCounts.txt

# counts from bowtie2 mapping (dual species)
srun featureCounts -T ${SLURM_CPUS_PER_TASK} -p -s 2 -O -M --fraction --primary -t gene -g ID -a mergeAnnot.gff -o all_bt2_featureCounts.txt *.hisat2_unmaped.bowtie2_sorted.bam
# srun head -n 2 all_bt2_featureCounts.txt | tail -n 1 | cut -d $'\t' -f 1,7-27 | sed 's/.hisat2_unmaped.bowtie2_sorted.bam//g' > all_bt2_mouse_featureCounts_matrix.tab
# srun head -n 2 all_bt2_featureCounts.txt | tail -n 1 | cut -d $'\t' -f 1,7-27 | sed 's/.hisat2_unmaped.bowtie2_sorted.bam//g' > all_bt2_cdiff_featureCounts_matrix.tab
# srun awk -F "\t" '{if(NR>2){if(($2!="NC_008226.2")&&($2!="NC_009089.1")){print $1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27>>"all_bt2_mouse_featureCounts_matrix.tab"}else{print $1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27>>"all_bt2_cdiff_featureCounts_matrix.tab"}}}' all_bt2_featureCounts.txt

echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
echo '########################################'
