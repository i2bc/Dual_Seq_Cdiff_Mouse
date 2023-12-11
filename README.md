# Dual_Seq_Cdiff_Mouse

Repository for the _C.difficiles_ RNAseq analyses of the Kreis _et al._ 2024 : "Dual RNA-seq study of the dynamics of coding and non-coding RNAs expression during Clostridioides difficile infection in a mouse model"

## Dual-RNAseq analysis

Dual-RNAseq data is available under identifier [PRJEB64651](https://www.ebi.ac.uk/ena/browser/view/PRJEB64651).

The RNAseq analysis process consists of several parts: quality control (01_initial_qc), cleaning (02_data_processing), mapping (03_mappin_genome), prediction of new candidate ncRNA genes (04_detrprok), counting the reads associated with each gene (05_feature_Counts_genome), analysis of differential gene expression (06_DESeq2), and comparisons (08_post_analysis_and_figures) with similar analyses already published (07_Fletcher_and_Pruss_analyses).

Install third-party softwares with the [conda](https://docs.conda.io/en/latest/) environment files present in the `00_initial_data/conda_env` repository (command line example: `conda env create -f 00_initial_data/conda_env/*.yml`)

### prediction of new candidate ncRNA genes

The prediction of new candidate ncRNA genes was done with [DETR'PROK, 2.1.3 version](http://rssf.i2bc.paris-saclay.fr/Software/detrprok.php):


### 07_Fletcher_and_Pruss_analyses

For this comparisons with similar published analyses part, the entire analysis workflow (01-06) was taken over in the form of a snakemake pipeline, as the data were not from a dual-RNAseq and therefore required adaptation. 

#### publisehd analysis used

- [Fletcher](https://doi.org/10.1038/s41467-020-20746-4), J.R., Pike, C.M., Parsons, R.J. et al. Clostridioides difficile exploits toxin-mediated inflammation to alter the host nutritional landscape and exclude competitors from the gut microbiota. Nat Commun 12, 462 (2021). https://doi.org/10.1038/s41467-020-20746-4
- [Pruss](https://doi.org/10.1038/s41586-021-03502-6), K.M., Sonnenburg, J.L. C. difficile exploits a host metabolite produced during toxin-mediated disease. Nature 593, 261–265 (2021). https://doi.org/10.1038/s41586-021-03502-6

List of sample ID and download url links stand in supplemental table 4.

#### get example data set

From the Pruss study, WT _vs._ IV conditions:
- WT: [SRR12762560](https://www.ebi.ac.uk/ena/browser/view/SRR12762560) [SRR12762561](https://www.ebi.ac.uk/ena/browser/view/SRR12762561) [SRR12762562](https://www.ebi.ac.uk/ena/browser/view/SRR12762562) [SRR12762563](https://www.ebi.ac.uk/ena/browser/view/SRR12762563)
- IV: [SRR12766943]https://www.ebi.ac.uk/ena/browser/view/SRR12766943 [SRR12766946]https://www.ebi.ac.uk/ena/browser/view/SRR12766946 [SRR12766947]https://www.ebi.ac.uk/ena/browser/view/SRR12766947
