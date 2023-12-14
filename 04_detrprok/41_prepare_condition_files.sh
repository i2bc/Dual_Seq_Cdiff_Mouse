#!/bin/bash
#
# Script to run DETR'PROK on Kreis (k), Prus (p), and Fletcher (f) data
# https://github.com/i2bc/Dual_Seq_Cdiff_Mouse
#
#
# Creation of condition files with the list of SRRs:
echo SRR11296256$'\n'SRR11296257$'\n'SRR11296259$'\n'SRR11296270$'\n'SRR11296281 > cond_fwtd2.txt
echo SRR11296276$'\n'SRR11296277$'\n'SRR11296278 > cond_fwtTY.txt
echo IVsampled_triplicate_merged  > cond_kIVsampled.txt
echo MI28H_triplicate_merged$'\n'MI32H_triplicate_merged > cond_kMI.txt
echo SRR12766943$'\n'SRR12766946$'\n'SRR12766947 > cond_pBase.txt
echo SRR12762560$'\n'SRR12762561$'\n'SRR12762562$'\n'SRR12762563 > cond_pWT.txt

