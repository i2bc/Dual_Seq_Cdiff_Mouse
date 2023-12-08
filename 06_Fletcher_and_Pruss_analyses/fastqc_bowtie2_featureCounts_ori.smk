# aim: add ncRNA to analyse the Fletcher RNAseq 2018 data
#
# may 2021
# usage example on IFB cluster (conda directives not tested):
# module load slurm-drmaa/1.0.8 snakemake/5.19.2 fastqc/0.11.9 samtools/1.10 bowtie2/2.4.1 subread/2.0.1
# sbatch snakemake --drmaa --jobs=4 -s fastqc_bowtie2_featureCounts.smk --configfile fastqc_bowtie2_featureCounts.yml 
# ############# ??? pour la création de l'index génomique bowtie2:
# sbatch --mem 80GB snakemake --drmaa --jobs=4 -s fastqc_bowtie2_featureCounts.smk --configfile fastqc_bowtie2_featureCounts.yml 

SAMPLES, = glob_wildcards(config["dataDir"]+"{sample}"+config["pePrefix"]+"1"+config["fastqSuff"])
BWT2_IDX = ["1","2","3","4","rev.1","rev.2"]
#CONDS, =   glob_wildcards("05_DEG/4Sugar/tables/{comparison}.complete.txt")

#totest:
#def getTargetFiles():
#    targets = list()
#    for r in config["st_dir"]:
#    targets.append([r]"data/"+config["refs"][r]+".fna.sa")
#    return targets


rule all:
  input:   
    #expand(config["fastqcDir"]+"{sample}"+config["pePrefix"]+"1_fastqc.html", sample=SAMPLES),
    #expand(config["fastqcDir"]+"{sample}"+config["pePrefix"]+"2_fastqc.html", sample=SAMPLES),
    #expand(config["bwt2Dir"]+config["genomeIdxPrefix"]+".{bwt2_idx}.bt2", bwt2_idx=BWT2_IDX),
    #expand(config["bwt2Dir"]+"{sample}.bam", sample=SAMPLES),
    #expand(config["bwt2Dir"]+"{sample}.flagstat", sample=SAMPLES),
    #expand(config["bwt2Dir"]+"{sample}.bam.bai", sample=SAMPLES),
    #expand(config["featureCountsDir"]+"{sample}_fc.txt", sample=SAMPLES),
    #expand(config["featureCountsDir"]+"{sample}_fc4DEG.txt", sample=SAMPLES),
    expand(config["cwdDir"]+config["st_dir"]+"{comparison}"+"/"+"{comparison}"+".RData", comparison=config["st_comparison"]),
    expand(config["cwdDir"]+config["st_dir"]+"{comparison}"+"/"+"{comparison}"+"_report.html", comparison=config["st_comparison"]),
    #
    expand(config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+".complete.txt"),
    expand(config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+".complete_annot.txt")
    #expand("05_DEG/4Sugar/tables/{condition}.complete_annot.txt", condition=config) 


rule addAnnotations:
  output:
    annot=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+".complete_annot.txt",
    tmp=temp(config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+".complete_tmp.txt")
    #annot="05_DEG/4Sugar/tables/{condition}.complete_annot.txt",
    #tmp=temp("05_DEG/4Sugar/tables/{condition}.complete_tmp.txt")
  input:
    #RData=config["st_dir"]+"{st_comparison}"+"/"+"{st_comparison}"+".RData"
    #dirAnalyseDiff=rules.analyseDiff.output.resultDir,
    #dirAnalyseDiff=config["st_dir"]+"{st_comparison}"+"/tables/",
    #table=config["st_dir"]+"{st_comparison}"+"/tables/"+"{st_condition}"+".complete.txt"
    table=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+".complete.txt"
  params:
    sep="-t $'\\t'",
    annot=config["st_designDir"]+config["annotFile"]
  shell:  
    """
    awk -f addAnnotationHeaderLine.awk {input.table} > {output.annot}
    sed 1d {input.table} | sort {params.sep} -k1,1 > {output.tmp}
    join {params.sep} {output.tmp} {params.annot} >> {output.annot} 
    """


rule analyseDiff:
  output: 
    RData=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/"+config["st_comparison"]+".RData",
    html=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/"+config["st_comparison"]+"_report.html",
#    resultDir=directory(config["st_dir"]+config["st_comparison"]+"/tables/")
    table=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+".complete.txt"
#,    table=config["st_dir"]+"{st_comparison}"+"/tables/"+lambda wildcards:config[wildcards.st_comparison]["st_conditions"]+".complete.txt"
  input: 
    design=config["cwdDir"]+config["st_designDir"]+config["st_comparison"]+".design4sartools",
    files=expand(config["featureCountsDir"]+"{sample}_fc4DEG.txt", sample=SAMPLES)
  params:
    resultDir=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/",
    html=config["st_comparison"]+"_report.html",
    RData=config["st_comparison"]+".RData",
    script=config["cwdDir"]+config["st_designDir"]+"sartools_script_DESeq2_CL.r",
    compName=config["st_comparison"],
    rawDir=config["cwdDir"]+config["featureCountsDir"],
    group=config["st_group"],
    condRef=config["st_condRef"]
  envmodules: "r/3.6.3", "pandoc/2.7.2"
  shell:
    """
    Rscript {params.script} --projectName="{params.compName}" --targetFile="{input.design}" --rawDir="{params.rawDir}" --varInt="{params.group}" --condRef="{params.condRef}"
    mv tables figures {params.html} {params.RData} {params.resultDir}
    """
   

rule prepareCounts:
  output: config["featureCountsDir"]+"{sample}_fc4DEG.txt"
  input: config["featureCountsDir"]+"{sample}_fc.txt"
  shell:
    """
    awk -f countFormating4DEG.awk {input} > {output}
    """

rule featureCount:
  output: config["featureCountsDir"]+"{sample}_fc.txt"
  input: config["bwt2Dir"]+"{sample}.bam"
  params:
    annotation=config["annotation"],
    feature=config["fc_feature"],
    tag=config["fc_tag"],
    mode=config["fc_mode"] 
  conda:
    config["condaEnvDir"]+"count.yml"
  envmodules: "subread/2.0.1"
  threads: 8
  shell:
    """
       srun featureCounts -T {threads} {params.mode} -t {params.feature} -g {params.tag} -a {params.annotation} -o {output} {input}
    """

rule bwt_bam_2_bai_and_stat:
  output:
    bai=config["bwt2Dir"]+"{sample}.bam.bai",
    fst=config["bwt2Dir"]+"{sample}.flagstat"
  input:
    config["bwt2Dir"]+"{sample}.bam"
  log:
    out_bai="Logs/{sample}_bwt_bam_2_bai.stdout",
    err_bai="Logs/{sample}_bwt_bam_2_bai.stderr",
    #out_stat="Logs/{sample}_bwt_bam_2_stat.stdout",
    err_stat="Logs/{sample}_bwt_bam_2_stat.stderr"
  conda:
    config["condaEnvDir"]+"samtools.yml"
  envmodules: "samtools/1.13"
  shell:
    """
       samtools index {input} {output.bai} 1> {log.out_bai} 2> {log.err_bai} ;
       samtools flagstat {input} 1> {output.fst} 2> {log.err_stat}
    """

rule bwt_sam_2_sorted_bam:
  output:
    config["bwt2Dir"]+"{sample}.bam"
  input:
    config["bwt2Dir"]+"{sample}.sam"
  log:
    out="Logs/{sample}_bwt_sam_2_sorted_bam.stdout",
    err="Logs/{sample}_bwt_sam_2_sorted_bam.stderr"
  conda:
    config["condaEnvDir"]+"samtools.yml"
  envmodules: "samtools/1.13"
  shell:
    "samtools sort -O BAM -o {output} {input} 1> {log.out} 2> {log.err}"


rule bowtie2_genome_indexing:
  output:
    expand(config["bwt2Dir"]+config["genomeIdxPrefix"]+".{ext}.bt2", ext=BWT2_IDX)
  input:
    fasta=config["genomeFna"]
  params:
    dir=config["bwt2Dir"]+config["genomeIdxPrefix"]
  resources:
#    cpus=8
    mem_mb=80000  
  threads : 8
  conda:
    config["condaEnvDir"]+"bowtie2.yml"
  envmodules: "bowtie2/2.4.1"
  log:
    out="Logs/genome_bwt2_index.stdout",
    err="Logs/genome_bwt2_index.stderr"
  shell:
    "bowtie2-build --threads {threads} {input.fasta} {params.dir} 1>{log.out} 2>{log.err}"


rule bowtie2_mapping:
  output:
    temp(config["bwt2Dir"]+"{sample}.sam")
  input:
    idx=rules.bowtie2_genome_indexing.output,
    r1=config["dataDir"]+"{sample}"+config["pePrefix"]+"1"+config["fastqSuff"],
    r2=config["dataDir"]+"{sample}"+config["pePrefix"]+"2"+config["fastqSuff"]
  params:
    index=config["bwt2Dir"]+config["genomeIdxPrefix"],
    orientation=config["bwtSequencingOrientation"]
#,    unmap=config["dataDir"]+"{sample}_bt2unmap.gz" # if option --un-conc-gz {params.unmap}
  log:
    out="Logs/{sample}_bwt2_mapping.stdout",
    err="Logs/{sample}_bwt2_mapping.stderr"
  resources: 
#    cpus=8,
    part="long",
    mem_mb=6000
  threads: 8
  conda: 
    config["condaEnvDir"]+"bowtie2.yml"
  envmodules: "bowtie2/2.4.1"
  shell:
    "bowtie2 -p {threads} --omit-sec-seq --sam-no-qname-trunc {params.orientation} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output} 1> {log.out} 2> {log.err}"


rule fastqc:
  output: 
    config["fastqcDir"]+"{sample}"+config["pePrefix"]+"1_fastqc.zip",
    config["fastqcDir"]+"{sample}"+config["pePrefix"]+"2_fastqc.zip",
    config["fastqcDir"]+"{sample}"+config["pePrefix"]+"1_fastqc.html",
    config["fastqcDir"]+"{sample}"+config["pePrefix"]+"2_fastqc.html"
  input: 
    r1=config["dataDir"]+"{sample}"+config["pePrefix"]+"1.fastq.gz",
    r2=config["dataDir"]+"{sample}"+config["pePrefix"]+"2.fastq.gz"
  params:
    outDir=config["fastqcDir"]
  log:
    out_r1="Logs/{sample}"+config["pePrefix"]+"1_fastqc.stdout",
    err_r1="Logs/{sample}"+config["pePrefix"]+"1_fastqc.stderr",
    out_r2="Logs/{sample}"+config["pePrefix"]+"2_fastqc.stdout",
    err_r2="Logs/{sample}"+config["pePrefix"]+"2_fastqc.stderr"
  conda: 
    config["condaEnvDir"]+"fastqc.yml"
  envmodules: "fastqc/0.11.9"
  shell: 
    """
       fastqc --outdir {params.outDir} {input.r1} 1>{log.out_r1} 2>{log.err_r1} ;
       fastqc --outdir {params.outDir} {input.r2} 1>{log.out_r2} 2>{log.err_r2}
    """
