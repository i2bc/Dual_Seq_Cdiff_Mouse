# Dual-RNASeq C.difficile & Mouse

Repository for the _C.difficiles_ RNAseq analyses of the Kreis _et al._ 2024 : "Dual RNA-seq study of the dynamics of coding and non-coding RNAs expression during Clostridioides difficile infection in a mouse model"

## Dual-RNAseq datd analysis

Dual-RNAseq data is available under identifier [PRJEB64651](https://www.ebi.ac.uk/ena/browser/view/PRJEB64651).

The RNAseq analysis process consists of several parts: quality control (`01_initial_qc`), cleaning (`02_data_preprocessing`), mapping (`03_mapping_genome`), prediction of new candidate ncRNA genes (`04_detrprok`), counting the reads associated with each gene (`05_featureCounts_genome`), analysis of differential gene expression (`06_DESeq2`), and comparisons (`08_post_analysis_and_figures`) with similar already published analyses (`07_Fletcher_and_Pruss_analyses`).

References genomes (fasta format) and annotations (gff format) used may be downloaded from:
- [Mus musculus genome assembly GRCm39](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27)
- [Clostridioides difficile 630](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009205.2). The annotation file for this organism was downloaded the 28th February 2023 from the [MicroScope](http://www.genoscope.cns.fr/agc/microscope) platform and given in the `00_initial_data` folder (`NC_009089.1_MaGe.gff`).

Prediction of new candidate ncRNA genes and comparisons with similar already published analyses concern only _C.difficiles_ organism.

The prediction of new candidate ncRNA genes was done with:
- [DETR'PROK v2.1.3](http://rssf.i2bc.paris-saclay.fr/Software/detrprok.php)
- `Clusterize v1.0.3` and `Compare Overlapping v1.0.4` of the [S-MART tools](https://urgi.versailles.inrae.fr/Tools/S-Mart)

### Fletcher and Pruss analyses

For the comparisons with similar published analyses part, the analysis workflow (01-06) was taken over in the form of a snakemake pipeline, as the data were not from a dual-RNAseq and therefore required adaptation. 

The publisehd analysis used are:

- [Fletcher](https://doi.org/10.1038/s41467-020-20746-4), J.R., Pike, C.M., Parsons, R.J. et al. Clostridioides difficile exploits toxin-mediated inflammation to alter the host nutritional landscape and exclude competitors from the gut microbiota. Nat Commun 12, 462 (2021). doi: 10.1038/s41467-020-20746-4
- [Pruss](https://doi.org/10.1038/s41586-021-03502-6), K.M., Sonnenburg, J.L. C. difficile exploits a host metabolite produced during toxin-mediated disease. Nature 593, 261–265 (2021). doi: 10.1038/s41586-021-03502-6

The files `selection_from_Fletcher_study.txt` and `selection_from_Pruss_study.txt` contain the selected samples used for the comparisons (see in `07_Fletcher_and_Pruss_analyses/data_example` folder).

Third-party softwares may be accessible with the [conda](https://docs.conda.io/en/latest/) environment files present in the `07_Fletcher_and_Pruss_analyses/conda_env` repository (command line example: `conda env create -f conda_env/*.yml`)

To test snakefile operation on a lightweight example:

1- go to the `07_Fletcher_and_Pruss_analyses` folder

2- download of the `template_script_DESeq2_CL.r` of the the [SARTools](https://github.com/PF2-pasteur-fr/SARTools) package in the `07_Fletcher_and_Pruss_analyses` folder (needed for the differential analysis step)

3- extract the archive `Dual_Seq_Cdiff_Mouse_smk_example.tar.gz` that contains the genome (fasta format) and annotation files (gff format) of _Clostridioides difficile 630_ and the first 10000 reads of R1 and R2 RNAseq data of the Pruss study (WT: [SRR12762560](https://www.ebi.ac.uk/ena/browser/view/SRR12762560) and [SRR12762561](https://www.ebi.ac.uk/ena/browser/view/SRR12762561) ; Base: [SRR12766943](https://www.ebi.ac.uk/ena/browser/view/SRR12766943) and [SRR12766946](https://www.ebi.ac.uk/ena/browser/view/SRR12766946)) in the `07_Fletcher_and_Pruss_analyses` folder: `tar -xvzf Dual_Seq_Cdiff_Mouse_smk_example.tar.gz`

4- adapt path of the `fQC_bwt2_ftCounts_DEseq2_annot.yml` config file to you own environment

5- Under the activated Dual_Seq_Cdiff_Mouse_smk conda environment (`conda activate env_Dual_Seq_Cdiff_Mouse_smk`) and in the `07_Fletcher_and_Pruss_analyses` folder, run : `snakemake --snakefile fQC_bwt2_ftCounts_DEseq2_annot.smk --configfile fQC_bwt2_ftCounts_DEseq2_annot.yml --cores 1`

6- The functional test is completed if there is no difference between the expected result file and the current result file (excluding columns based on statistical estimation): `diff <(cut -f 1-12,21- data_example/expectedResultOfFunctionalTest.txt) <(cut -f 1-12,21- 05_DEG/functional_example/tables/WTvsBase.complete_annot.txt)` (`07_Fletcher_and_Pruss_analyses` folder)

### Post Analysis and Figures

Post analysis and figures for the three experiments are based on the tables from the differential expression analyses for all genes (`complete`) for each experiment. Dedicated files and scripts stand in the `08_post_analysis_and_figures` directory. 

Chi2 test comparisons and heapmaps (with the [ComplexHeatmap](http://bioconductor.org/packages/ComplexHeatmap/) package) of differentially expressed genes in at least two experiments are performed by the `figures_stats_noNA4allStat.R` script. The corresponding [conda](https://docs.conda.io/en/latest/) environment may be created with the `conda_env/env_R_for_figures_and_stats_Kreis.yml` file.

Gene set enrichment analyses (GSEA) were performed using a manual replay (`Ma2Html_for_blitzgsea_analysis`) of the [Ma2HTML](https://mmonot.eu/MA2HTML/connexion.php)(M. Monot personnal communication, export 1652263129) database to define gene sets, the Wald statistic (`stat` column proposed by [SARTools](https://github.com/PF2-pasteur-fr/SARTools)) as gene signatures, and the [blitzgsea](https://github.com/MaayanLab/blitzgsea/) tool using the `blitzgsea_Ma2Html.py` script. 

The pie charts are based on the differentially expressed genes in the leading_edge lists from the GSEA analysis.



