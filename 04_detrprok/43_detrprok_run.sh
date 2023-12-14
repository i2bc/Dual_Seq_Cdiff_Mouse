#!/bin/bash
#
# Script to run DETR'PROK on Kreis (k), Prus (p), and Fletcher (f) data
# https://github.com/i2bc/Dual_Seq_Cdiff_Mouse
#
# Thrid-party softwares:
# - DETRPROK_2.1.3.sh http://rssf.i2bc.paris-saclay.fr/Software/detrprok.php
#
# DETR'PROK runs with 2 sets of parameters following experiments
mkdir detrprok_res ; 
covaS=100 ; covsR=3 ; cov5u=5 ; for i in `cat cond_fwtd2.txt cond_fwtTYsampled.txt cond_pBaseSampled cond_pWT.txt ` ; do 
   bash DETRPROK_2.1.3.sh -bed bam_bed/${i}_ppR.bed -gff ../00_initial_data/NC_009089.1_dnaA.gff -read_len 100 -features_list "CDS|ncRNA|riboswitch|RNase_P_RNA|rRNA|SRP_RNA|tmRNA|tRNA" -op_gap 60 -clust_gap 20 -RNA_gap 25 -RNA_merge 25 -sRNA_coverage ${covsR} -asRNA_coverage ${covaS} -5utr_coverage ${cov5u} -verbose true ;
   for j in 5UTRs sRNAs asRNAs ; do 
      mv bam_bed/${i}_ppR_${j}.gff detrprok_res/${i}_ppR_${j}_opgap60_clustgap20_RNAgap25_RNAmerge25_covsR${covsR}_cov5u${cov5u}_covaS${covaS}.gff ; 
   done ; 
done 
covaS=10 ; covsR=2 ; cov5u=5 ; for i in `cat cond_kIVsampled.txt cond_kMI.txt ` ; do 
   bash DETRPROK_2.1.3.sh -bed bam_bed/${i}_ppR.bed -gff ../00_initial_data/NC_009089.1_dnaA.gff -read_len 100 -features_list "CDS|ncRNA|riboswitch|RNase_P_RNA|rRNA|SRP_RNA|tmRNA|tRNA" -op_gap 60 -clust_gap 20 -RNA_gap 25 -RNA_merge 25 -sRNA_coverage ${covsR} -asRNA_coverage ${covaS} -5utr_coverage ${cov5u} -verbose true ;
   for j in 5UTRs sRNAs asRNAs ; do 
      mv bam_bed/${i}_ppR_${j}.gff detrprok_res/${i}_ppR_${j}_opgap60_clustgap20_RNAgap25_RNAmerge25_covsR${covsR}_cov5u${cov5u}_covaS${covaS}.gff ; 
   done ; 
done 
