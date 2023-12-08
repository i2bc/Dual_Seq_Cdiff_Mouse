#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=631_adding_mouseGeneID

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

# usage example:
# sbatch 63_adding_product.sh early_late_cdiff/tables/earlyvslate.complete.txt
# limits:
# - not for "up" or "down" DESeq2 tables

# Main script
count_file=`basename $1 .txt`
count_dir=`dirname $1`

# determine species of interest:
#    it is cdiff species if grep -c > 0
   # Mus species
   annot=mouse_geneID.list
   if [ ! -f ${annot} ]; then 
      # create annotations list of 2 columns (gene name from Parent field + GeneID included in the Dbxref field)
      mus_annots=/shared/projects/clos_ribo_reg/reference_dbs/Mmusculus/GRCm39/GCF_000001635.27_GRCm39_genomic.gff
      srun awk -F "\t" -v id="ID=" -v dbxref="Dbxref=" -v geneID="GeneID" '{if($3=="gene"){idVal[2]="";dbxrefVal[2]="";nb=split($9,col9,";");for(f=1;f<=nb;f++){if(col9[f]~id){split(col9[f],idVal,"=")};if(col9[f]~dbxref){split(col9[f],dbxrefVal,"=");nb_ref=split(dbxrefVal[2],ref,",");for(x=1;x<=nb_ref;x++){if(ref[x]~geneID){split(ref[x],ID,":")}}}};printf("%s\t%s\n",idVal[2],ID[2])}}' ${mus_annots} | sort -t $'\t' -k1,1 > mouse_geneID.list
   fi

# add annotations to the SARtools "complete" table
# i) header line:
awk 'NR==1{print(sprintf("%s\tGeneID",$0))}' ${count_dir}/${count_file}.txt > ${count_dir}/${count_file}_annot_4david.tsv
# ii) prepare counts table (suppress header + sort):
sed 1d ${count_dir}/${count_file}.txt | sort -t $'\t' -k1,1 > ${count_dir}/${count_file}_tmp.txt
# iii) add annotations to counts tables:
join -t $'\t' ${count_dir}/${count_file}_tmp.txt ${annot} >> ${count_dir}/${count_file}_annot_4david.tsv
rm ${count_dir}/${count_file}_tmp.txt

# select columns (Id "norm.x" cond1 cond2 log2FoldChange padj GeneID4david):
#    minMorm = 1rst column with norm.x
srun awk -F "\t" 'NR==1{minNorm=NF;for(c=2;c<=NF-14;c++){if((match($c,"norm."))&&(minNorm>c)){minNorm=c}}}NR>=1{printf("%s",$1);for(i=minNorm;i<=NF-15;i++){printf("\t%s",$i)};printf("\t%s\t%s\t%s\t%s\t%s\n",$(NF-13),$(NF-12),$(NF-10),$(NF-7),$NF)}' ${count_dir}/${count_file}_annot_4david.tsv > ${count_dir}/${count_file}_annot_4david_reduced.tsv
