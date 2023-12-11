#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=61_fcMatrix_2_SamplesCounts

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

# Main script
dir_fcMatrix=../05_featureCounts_genome
rm fc_data/all_hbt2_featureCounts_matrix.tab
rm fc_data/*_CD630_fc.txt 
rm fc_data/*_mouse_fc.txt
# create matrix with ht2 + bt2 :
  # header line:
awk -F "\t" 'NR==2{printf("%s",$1);for(i=7;i<=NF;i++){printf("\t%s",$i)};printf("%s","\n")}' ${dir_fcMatrix}/all_ht2_featureCounts.txt | sed 's/.cleaned.hisat2_sorted.bam//g' > fc_data/all_hbt2_featureCounts_matrix.tab ; 
  # count lines:
paste ${dir_fcMatrix}/all_bt2_featureCounts.txt ${dir_fcMatrix}/all_ht2_featureCounts.txt | awk -F "\t" 'NR>2{printf("%s",$1);for(i=7;i<=(NF/2);i++){printf("\t%d",int($i+$(i+26)))};printf("%s","\n")}' >> fc_data/all_hbt2_featureCounts_matrix.tab
# create all species files :
  # tRNA suppression: mouse => 395 gene-n-T + 27 gene-Trn ; cdiff => 87 CD630_t
  # rRNA suppression: mouse => 51 gene-n-R + 2 gene-Rn + 2 gene-mt-R + 9 others ; cdiff => 87 CD630_t
  # how to find gene list : 
  # awk -F "\t" '{if($3=="gene"){print $0}}' GCF_000001635.27_GRCm39_genomic.gff | grep "gene_biotype=rRNA" | grep -v gene-n-R | grep -v gene-Rn | grep -v gene-mt-R | wc -l
#GeneToExclude="gene-n-T,gene-Trn,CD630_t,CD630_r,gene-n-R,gene-mt-R,gene-Rn4.5s,gene-Rn5s,gene-LOC115490001,gene-Gm25018,gene-Gm24312,gene-Gm26391,gene-Gm22109,gene-Gm25212,gene-Gm22291,gene-Gm23284,gene-LOC115488082"

# head fc_data/all_hbt2_featureCounts_matrix.tab
# Geneid	IV1_S3	IV1_S3_sampled	IV2_S4	IV2_S4_sampled	IV3_S5	IV3_S5_sampled	S10_S5	S11_S1	S12_S2	S1_S2	S2_S3	S3_S4	S4_S5	S5_S6	S6_S1	S7_S2	S8_S3	S9_S46	S9_S4	S9_S6

awk -F "\t" -v exclude=$(cat genesToExclude.txt) 'BEGIN{nbex=split(exclude,ex,",")}NR==1{for(i=2;i<=NF;i++){Sa[i]=$i}};NR>1{exclusion=0;for(e=1;e<=nbex;e++){if($1==ex[e]){exclusion=1}};if(exclusion==0){if($1~/CD630_/){species="cdiff"}else{species="mouse"};for(i=2;i<=NF;i++){print $1"\t"$i>>"fc_data/"Sa[i]"_"species"_fc.txt"}}}' fc_data/all_hbt2_featureCounts_matrix.tab
echo '########################################'
echo 'Job finished' $(date --iso-8601=seconds)
