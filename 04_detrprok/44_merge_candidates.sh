#!/bin/bash
#
# Script to run DETR'PROK on Kreis (k), Prush (p), and Fletcher (f) data
# https://github.com/i2bc/Dual_Seq_Cdiff_Mouse
#
# Thrid-party software:
# - clusterize.py and CompareOverlapping.py of the S-MART tool (python 2.7 writen) https://urgi.versailles.inrae.fr/Tools/S-Mart 
#
# Associate candidates
rm sRNA_candidates.tmp ; 
for i in kMI kIVsampled pWT pBaseSampled fwtTYsampled fwtd2 ; do 
   for s in `cat cond_${i}.txt` ; do 
      sed "s/sRNA_/${i}_sR_/" detrprok_res/${s}_ppR_sRNAs_*.gff >> sRNA_candidates.tmp ; 
      sed "s/asRNA_/${i}_asR_/" detrprok_res/${s}_ppR_asRNAs_*.gff >> sRNA_candidates.tmp ; 
   done
done
Smart/clusterize.py -i sRNA_candidates.tmp -f gff -o sRNA_candidates.tmp1 -u gff -c 
Smart/CompareOverlapping.py -j ../00_initial_data/sRNA_knownBeforeThisWork.gff -f gff -i sRNA_candidates.tmp1 -g gff -c -x -o sRNA_candidates.tmp2
Smart/CompareOverlapping.py -j ../00_initial_data/NC_009089.1_MaGe_rtRNA.gff -f gff -i sRNA_candidates.tmp2 -g gff -d 100 -c -x -o sRNA_candidates.tmp3
Smart/CompareOverlapping.py -j ../00_initial_data/NC_009089.1_MaGe_rtRNA.gff -f gff -i sRNA_candidates.tmp3 -g gff -a -x -o sRNA_candidates.tmp4
sed 's/Name=/detectedIn=/;s/S-MART/detrprok/;s/asRNA/sRNA/' sRNA_candidates.tmp4 > sRNA_candidates.tmp5
deleteTagGff.pl -i sRNA_candidates.tmp5 -d nbElements,nbOverlappingReads,ID > sRNA_candidates.gff

