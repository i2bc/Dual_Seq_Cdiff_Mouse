# Dual-RNASeq C.difficile & Mouse

Repository for the _C.difficiles_ RNAseq analyses of the Kreis _et al._ 2024 : "Dual RNA-seq study of the dynamics of coding and non-coding RNAs expression during Clostridioides difficile infection in a mouse model"

## Dual-RNAseq datd analysis

Dual-RNAseq data is available under identifier [PRJEB64651](https://www.ebi.ac.uk/ena/browser/view/PRJEB64651).

The RNAseq analysis process consists of several parts: quality control (`01_initial_qc`), cleaning (`02_data_preprocessing`), mapping (`03_mapping_genome`), prediction of new candidate ncRNA genes (`04_detrprok`), counting the reads associated with each gene (`05_featureCounts_genome`), analysis of differential gene expression (`06_DESeq2`), and comparisons (`08_post_analysis_and_figures`) with similar already published analyses (`07_Fletcher_and_Pruss_analyses`).

References genomes (fasta format) and annotations (gff format) used may be downloaded from:
- [Mus musculus genome assembly GRCm39](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27)
- [Clostridioides difficile 630](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009205.2)

Prediction of new candidate ncRNA genes and comparisons with similar already published analyses concern only _C.difficiles_ organism.

### Prediction of new candidate ncRNA genes

The prediction of new candidate ncRNA genes was done with:
- [DETR'PROK v2.1.3](http://rssf.i2bc.paris-saclay.fr/Software/detrprok.php)
- `Clusterize v1.0.3` and `Compare Overlapping v1.0.4` of the [S-MART tools](https://urgi.versailles.inrae.fr/Tools/S-Mart)


### Fletcher and Pruss analyses

For the comparisons with similar published analyses part, the analysis workflow (01-06) was taken over in the form of a snakemake pipeline, as the data were not from a dual-RNAseq and therefore required adaptation. 

The publisehd analysis used are:

- [Fletcher](https://doi.org/10.1038/s41467-020-20746-4), J.R., Pike, C.M., Parsons, R.J. et al. Clostridioides difficile exploits toxin-mediated inflammation to alter the host nutritional landscape and exclude competitors from the gut microbiota. Nat Commun 12, 462 (2021). doi: 10.1038/s41467-020-20746-4
- [Pruss](https://doi.org/10.1038/s41586-021-03502-6), K.M., Sonnenburg, J.L. C. difficile exploits a host metabolite produced during toxin-mediated disease. Nature 593, 261–265 (2021). doi: 10.1038/s41586-021-03502-6

The files `selection_from_Fletcher_study.txt` and `selection_from_Pruss_study.txt` contain the selected samples used for the published comparisons (see in `07_Fletcher_and_Pruss_analyses/data_example` folder).

Third-party softwares may be accessible with the [conda](https://docs.conda.io/en/latest/) environment files present in the `00_initial_data/conda_env` repository (command line example: `conda env create -f 00_initial_data/conda_env/*.yml`)

To run the snakemake pipeline on a toy example :

The genome (fasta format) and annotation files (gff format) of _Clostridioides difficile 630_ should be present in the `07_Fletcher_and_Pruss_analyses/data_example` repository.

The differential analysis step is done with the [SARTools](https://github.com/PF2-pasteur-fr/SARTools) package. Copy of the `template_script_DESeq2_CL.r` in the `07_Fletcher_and_Pruss_analyses` folder.

From the Pruss study, WT _vs._ IV conditions, download R1 and R2 RNAseq data:
- WT: [SRR12762560](https://www.ebi.ac.uk/ena/browser/view/SRR12762560) [SRR12762561](https://www.ebi.ac.uk/ena/browser/view/SRR12762561)
- IV: [SRR12766943](https://www.ebi.ac.uk/ena/browser/view/SRR12766943) [SRR12766946](https://www.ebi.ac.uk/ena/browser/view/SRR12766946) 







