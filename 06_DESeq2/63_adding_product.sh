#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=63_adding_product

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
test_species=`grep -c CD630 ${count_dir}/${count_file}.txt` 
if [ ${test_species} -gt 0 ]; then
   # Cdiff species :
   # create or fix anotation file
   annot=cdiff_id_prod.list
   if [ ! -f ${annot} ]; then # -f => existing file?
      # create annotations list of 2 columns (gene name from Parent field + product)
      cdiff_annots=/shared/projects/clos_ribo_reg/reference_dbs/Cdifficile_MaGe/CD630_MaGe_ID.gff
      srun awk -F "\t" -v id=Parent -v prod=product '{if(($3=="CDS")||($3=="misc_RNA")){idVal[2]="";prodVal[2]="";nb=split($9,col9,";");for(i=1;i<=nb;i++){if(col9[i]~id){split(col9[i],idVal,"=")};if(col9[i]~prod){split(col9[i],prodVal,"=")}};printf("%s\t%s\n",idVal[2],prodVal[2])}}' ${cdiff_annots} | sort -t $'\t' -k1,1 > cdiff_id_prod.list
   fi
else
   # Mus species
   annot=mouse_id_prod.list
   if [ ! -f ${annot} ]; then 
      # create annotations list of 2 columns (gene name from Parent field + product from description or biotype field)
      mus_annots=/shared/projects/clos_ribo_reg/reference_dbs/Mmusculus/GRCm39/GCF_000001635.27_GRCm39_genomic.gff
      srun awk -F "\t" '{if($3=="gene"){print}}' ${mus_annots} | awk -F "\t" -v id="ID=" -v prod="description=" -v type="gene_biotype=" '{idVal[2]="";prodVal[2]="";typeVal[2]="";nb=split($9,col9,";");for(i=1;i<=nb;i++){if(col9[i]~id){split(col9[i],idVal,"=")};if(col9[i]~prod){split(col9[i],prodVal,"=")};if(col9[i]~type){split(col9[i],typeVal,"=")}};if(prodVal[2]==""){prodVal[2]=typeVal[2]};printf("%s\t%s\n",idVal[2],prodVal[2])}' | sort -t $'\t' -k1,1 > mouse_id_prod.list
   fi
fi

# add annotations to the SARtools "complete" table
# i) header line:
awk 'NR==1{print(sprintf("%s\tproduct",$0))}' ${count_dir}/${count_file}.txt > ${count_dir}/${count_file}_annot.tsv
# ii) prepare counts table (suppress header + sort):
sed 1d ${count_dir}/${count_file}.txt | sort -t $'\t' -k1,1 > ${count_dir}/${count_file}_tmp.txt
# iii) add annotations to counts tables:
join -t $'\t' ${count_dir}/${count_file}_tmp.txt ${annot} >> ${count_dir}/${count_file}_annot.tsv
rm ${count_dir}/${count_file}_tmp.txt

# select columns (Id "norm.x" cond1 cond2 log2FoldChange padj product):
#    minMorm = 1rst column with norm.x
srun awk -F "\t" 'NR==1{minNorm=NF;for(c=2;c<=NF-14;c++){if((match($c,"norm."))&&(minNorm>c)){minNorm=c}}}NR>=1{printf("%s",$1);for(i=minNorm;i<=NF-15;i++){printf("\t%s",$i)};printf("\t%s\t%s\t%s\t%s\t%s\n",$(NF-13),$(NF-12),$(NF-10),$(NF-7),$NF)}' ${count_dir}/${count_file}_annot.tsv > ${count_dir}/${count_file}_annot_reduced.tsv
