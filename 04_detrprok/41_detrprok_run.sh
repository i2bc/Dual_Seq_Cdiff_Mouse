#!/bin/bash
#
# Script to run DETR'PROK on Kreis (k), Prus (p), and Fletcher (f) data
# https://github.com/i2bc/Dual_Seq_Cdiff_Mouse
#
#
cd 04_detrprok
# Creation of condition files with the list of SRRs:
echo SRR11296256$'\n'SRR11296257$'\n'SRR11296259$'\n'SRR11296270$'\n'SRR11296281 > cond_fwtd2.txt
echo SRR11296276$'\n'SRR11296277$'\n'SRR11296278 > cond_fwtTY.txt
echo IVsampled_triplicate_merged  > cond_kIVsampled.txt
echo MI28H_triplicate_merged$'\n'MI32H_triplicate_merged > cond_kMI.txt
echo SRR12766943$'\n'SRR12766946$'\n'SRR12766947 > cond_pBase.txt
echo SRR12762560$'\n'SRR12762561$'\n'SRR12762562$'\n'SRR12762563 > cond_pWT.txt
#
# Manage paired-end reads for DETR'PROK analysis:
# - selection of forward (reverse) primary mapping on chr NC_009089.1 from bam (samtools)
# - sampling of too large bam (for fwtTY and pBase conditions)
# - switch to bed format (bamToBed)
# - inverse the strand of the R1 file (awk)
# - concatenate the reversed R1 file and the R2 file (cat)
# - sort by ascending position for bedtools (sort)
conda activate samtools
mkdir bam_bed  ; 
for i in 03_mapping_genome/*.bam ; do samtools index ${i} ; done
for s in `cat cond_*.txt ` ; do 
   samtools view -f 131 -h -o bam_bed/${s}_ppR2.bam 03_mapping_genome/${s}.bam NC_009089.1 ; 
   samtools view -f  67 -h -o bam_bed/${s}_ppR1.bam 03_mapping_genome/${s}.bam NC_009089.1 ; 
done
rm cond_fwtTYsampled.txt ;
for s in `cat cond_fwtTY.txt ` ; do 
   samtools view -f 131 -s 11.2 -h -o bam_bed/${s}_sampled_ppR2.bam bam_bed/${s}.bam NC_009089.1 ; 
   samtools view -f  67 -s 11.2 -h -o bam_bed/${s}_sampled_ppR1.bam bam_bed/${s}.bam NC_009089.1 ; 
   cat ${s}_sampled >> cond_fwtTYsampled.txt ;
done
rm cond_pBaseSampled.txt ; 
for s in `cat cond_pBase.txt ` ; do 
   samtools view -f 131 -s 12.1 -h -o bam_bed/${s}_sampled_ppR2.bam bam_bed/${s}.bam NC_009089.1 ; 
   samtools view -f  67 -s 12.1 -h -o bam_bed/${s}_sampled_ppR1.bam bam_bed/${s}.bam NC_009089.1 ; 
      cat ${s}_sampled >> cond_pBaseSampled.txt ;
done
conda deactivate ; conda activate bedtools
for s in `cat cond_fwtd2.txt cond_fwtTYsampled.txt cond_kIVsampled.txt cond_kMI.txt cond_pBaseSampled.txt cond_pWT.txt ` ; do 
   for p in R1 R2 ; do 
         bamToBed -i bam_bed/${s}_pp${p}.bam > bam_bed/${s}_pp${p}.bed
   done
   awk 'BEGIN{FS="\t";OFS="\t"}{if($6=="-"){strand="+"}else{strand="-"};print $1,$2,$3,$4,$5,strand}' bam_bed/${s}_ppR1.bed > bam_bed/${s}_ppR.tmp 
   cat < bam_bed/${s}_ppR2.bed >> bam_bed/${s}_ppR.tmp
   sort -t $'\t' -k1,1 -nk2,2 -nk3,3 bam_bed/${s}_ppR.tmp > bam_bed/${s}_ppR.bed
done 
conda deactivate
#
# 2 DETR'PROK runs with 2 sets of parameters 
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
#
# Associate candidates
conda activate python2.7_ce
rm cand_tmp.gff ; 
for i in kMI kIVsampled ; do sed "s/sRNA_/${i}_sR_/" detrprok_res/${i}_sRNA.gff >> cand_tmp.gff ; sed "s/asRNA_/${i}_asR_/" detrprok_res/${i}*_asRNAs_*_covsR2_cov5u3_covaS10.gff >> cand_tmp.gff ; done
for i in pWT pBaseSampled fwtTYsampled fwtd2 ; do sed "s/sRNA_/${i}_sR_/" detrprok_res/${i}_sRNA.gff >> cand_tmp.gff ; sed "s/asRNA_/${i}_asR_/" detrprok_res/${i}*_asRNAs_*_covsR3_cov5u5_covaS100.gff >> cand_tmp.gff ; done
python ~/Sources/Smart/clusterize.py -i cand_tmp.gff -f gff -o cand_tmp_clust.gff -u gff -c 
grep nbE cand_tmp_clust.gff > cand_tmp_clust_only.gff
cp cand_tmp_clust_only.gff sRNA_candidates.tmp1
awk -F "\t" '{if(($3=="sRNA")&&($0!~"nbE")){print}}' cand_tmp_clust.gff >> sRNA_candidates.tmp1
Smart/CompareOverlapping.py -j ../00_initial_data/sRNA_knownBeforeThisWork.gff -f gff -i sRNA_candidates.tmp1 -g gff -c -x -o sRNA_candidates.tmp2
Smart/CompareOverlapping.py -j NC_009089.1_stRNA.gff -f gff -i sRNA_candidates.tmp2 -g gff -d 100 -c -x -o sRNA_candidates.tmp3
Smart/CompareOverlapping.py -j NC_009089.1_stRNA.gff -f gff -i sRNA_candidates.tmp3 -g gff -a -x -o sRNA_candidates.tmp4
deleteTagGff.pl -i sRNA_candidates.tmp4  -d nbElements,nbOverlappingReads,Name,ID | sed 's/Name=/detectedIn=/;s/S-MART/detrprok/;s/asRNA/sRNA/' > sRNA_candidates.gff
cd ..

