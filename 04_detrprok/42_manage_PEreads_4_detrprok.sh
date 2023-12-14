#!/bin/bash
#
# Script to run DETR'PROK on Kreis (k), Prus (p), and Fletcher (f) data
# https://github.com/i2bc/Dual_Seq_Cdiff_Mouse
#
# Thrid-party softwares:
# - samtools https://github.com/samtools/samtools 
# - bamToBed from bedtools https://bedtools.readthedocs.io/
#
# Manage paired-end reads for DETR'PROK analysis:
# - selection of forward (reverse) primary mapping on chr NC_009089.1 from bam (samtools)
# - sampling of too large bam (for fwtTY and pBase conditions)
# - switch to bed format (bamToBed)
# - inverse the strand of the R1 file (awk)
# - concatenate the reversed R1 file and the R2 file (cat)
# - sort by ascending position for next bedtools step (sort)
for i in ../03_mapping_genome/*.bam ; do samtools index ${i} ; done
mkdir bam_bed  ; 
for s in `cat cond_*.txt ` ; do 
   samtools view -f 131 -h -o bam_bed/${s}_ppR2.bam ../03_mapping_genome/${s}.bam NC_009089.1 ; 
   samtools view -f  67 -h -o bam_bed/${s}_ppR1.bam ../03_mapping_genome/${s}.bam NC_009089.1 ; 
   for r in 1 2 ; do samtools index bam_bed/${s}_ppR${r}.bam ; done      
done
rm cond_fwtTYsampled.txt ;
for s in `cat cond_fwtTY.txt ` ; do 
   samtools view -f 131 -s 11.2 -h -o bam_bed/${s}_sampled_ppR2.bam bam_bed/${s}_ppR2.bam NC_009089.1 ; 
   samtools view -f  67 -s 11.2 -h -o bam_bed/${s}_sampled_ppR1.bam bam_bed/${s}_ppR1.bam NC_009089.1 ; 
   echo ${s}_sampled >> cond_fwtTYsampled.txt ;
done
rm cond_pBaseSampled.txt ; 
for s in `cat cond_pBase.txt ` ; do 
   samtools view -f 131 -s 12.1 -h -o bam_bed/${s}_sampled_ppR2.bam bam_bed/${s}_ppR2.bam NC_009089.1 ; 
   samtools view -f  67 -s 12.1 -h -o bam_bed/${s}_sampled_ppR1.bam bam_bed/${s}_ppR1.bam NC_009089.1 ; 
   echo ${s}_sampled >> cond_pBaseSampled.txt ;
done
for s in `cat cond_fwtd2.txt cond_fwtTYsampled.txt cond_kIVsampled.txt cond_kMI.txt cond_pBaseSampled.txt cond_pWT.txt ` ; do 
   for p in R1 R2 ; do 
         bamToBed -i bam_bed/${s}_pp${p}.bam > bam_bed/${s}_pp${p}.bed
   done
   awk 'BEGIN{FS="\t";OFS="\t"}{if($6=="-"){strand="+"}else{strand="-"};print $1,$2,$3,$4,$5,strand}' bam_bed/${s}_ppR1.bed > bam_bed/${s}_ppR.tmp 
   cat < bam_bed/${s}_ppR2.bed >> bam_bed/${s}_ppR.tmp
   sort -t $'\t' -k1,1 -nk2,2 -nk3,3 bam_bed/${s}_ppR.tmp > bam_bed/${s}_ppR.bed
done 

