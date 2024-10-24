## RNA-Seq Pipeline and Analysis

<details open>
  <summary><b>Introduction</b></summary>


>  RNA-Seq is a sequencing method used to determine gene expression levels. [RNA-Seq](https://pmc.ncbi.nlm.nih.gov/articles/PMC6096346/) data originates from extracted RNA that was reverse transcribed into DNA.
There are many readyly available RNA-seq pipeline workflow eg. [nf-core](https://nf-co.re/rnaseq/3.14.0/), [PiGx](https://bioinformatics.mdc-berlin.de/pigx_docs/pigx-rna-seq.html), etc. [Galaxy](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html) is a popular framework software for constructing and managing the execution of workflows. 
  
</details>

<details>
  <summary><b>1 Tools & Workflow </b></summary>
  
#### 1.1 List of Tools:

>  List of tools commonly used in each step of RNA-seq data analysis. Depending on the the requirements select the tools accordingly. 

1. Quality Control (QC): [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://github.com/MultiQC/MultiQC), & [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
2. Improving Quality (Trimming and Filtering): [fastp](https://github.com/OpenGene/fastp), [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic), [Cutadapt](https://cutadapt.readthedocs.io/en/stable/), & [Skewer](https://github.com/relipmoc/skewer)
3. Read Alignment: [HISAT2](https://daehwankimlab.github.io/hisat2/manual/), [STAR](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html), [TopHat2](https://ccb.jhu.edu/software/tophat/manual.shtml), [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
4. Transcript Assembly and Quantification: [StringTie](), [Cufflinks](), [Salmon](), [Kallisto]()
5. Differential Expression Analysis: [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), enrichplot, pathview, edgeR, limma- R packages
6. Gene Ontology (GO) and KEGG - Pathway Analysis: [GSEA (Gene Set Enrichment Analysis)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page), clusterProfiler, EnhancedVolcano, DAVID
7. Visualization: [IGV (Integrative Genomics Viewer)](https://igv.org/doc/desktop/), pheatmap, ggplot2(R packages)
8. Single-Cell RNA-seq Specific Tools: CellRanger, Seurat, & SC3

***

#### 1.2 Workflow used for this Analysis :

>  For this RNA-seq analysis, **_FastQC_** is used for initial quality control, with **_MultiQC_** summarizing the results. **_fastp_** trims adapters and low-quality bases to improve read quality before alignment with **_HISAT2_**. Differential expression analysis is done using **_DESeq2_** R package. For gene ontology (GO) and pathway analysis (KEGG), **_GSEA_** and **_clusterProfiler_** are employed, with **_pathview_** used to visualize pathways. Visualization tools include **_enrichplot, emapplot,EnhancedVolcano_** for enrichment and DEG results, while **_pheatmap, ggplot2_** handle heatmaps and other graphical representations of the data.


***

</details>
  
<details>
  <summary><b>2. Installation Guide</b></summary>

  
  2.1 System Info:
  - **System:** _Ubuntu 24.04 LTS_ `lsb_release -a`
  - **RAM - threads:** _7.45G - 11threads_  `htop`
  - **Specs:** 172G avail `df -h`

  ---
  
  2.2 Conda Env Set-up:

- Install [Miniconda](https://docs.anaconda.com/miniconda/)
- Add [Bioconda channels](https://bioconda.github.io/)
- Create env:
```
conda create -n ranaseq
conda activate rnaseq
```
- Installation:
  - QC Tools:
  - Alignment Tools:
```
conda install bioconda::multiqc
conda install bioconda:fastqc
conda install -c bioconda fastp
```
 
```
conda install bioconda::samtools
conda install bioconda::hisat2
conda install bioconda::subread
```

---
     
  2.3 Conda Env: once installed activate the env `conda activate rnaseq_env`

```
conda env create -f rnaseq_env.yml
```

---

  2.4 Tools versions: use tool_name `--version` or  `-v` to the version `--help` or `-h` for user guide of the tool.
     
|Sr.no|Tools|Version|
|:----|:----|:-----:|
|1. |Miniconda|conda 24.5.0|
|2. |Python|Python 3.12.4|
|3. |FastQC|FastQC v0.12.1|
|4. |MultiQC|version 1.18|
|5. |HISAT2|version 2.2.1|
|6. |samtools|samtools 1.19.2|
|7. |Subread|featureCounts v2.0.1|

---
  </details>


<details>
  <summary><b>3 Pipeline Scripts & Automation </b></summary>
  
  3.1 Run commands: give flg links
  - build genome
  -
  -
  -
---
  
  3.2 Scripts:
  - QC Script:
  - Alignment Script:

---
  
  3.3 R script

---

  3.3 Automation:
  
- bash: things to add checkers if file present or not if present skip and start nex step
- lines to change eg in counts.sh change script/qc.sh and hisat2.sh path before running and nakesure your genome is built with the name genome if not change it xyz name you have ion line no.
```
#run "bash path_to_script_folder/count.sh path_to_rawfastq_folder/ path_to_genome_folder/ path_to_gtf-gff_file/.gtf
# main output of this script is featurecounts.tsv
bash scripts/count.sh ../biostateai/raw_fastq ../mgiGenome ../mgiGenome/gencode.vM35.basic.annotation.gtf 
```
- Folder Structure:

  biostateai/
├── scripts
│   ├── count.sh
│   ├── hisat2.sh
│   └── qc.sh
└── raw_fastq
    ├── bam
    │   ├── all_bam.txt
    │   ├── all_bam.txt.summary
    |   ├── featurecounts.tsv
    │   ├── Liver_ZT0_1.bam
    │   └── Liver_ZT12_1.bam
    ├── fastp
    │   ├── Liver_ZT0_1_fastp_error.log
    │   ├── Liver_ZT0_1_fastp.html
    │   ├── Liver_ZT0_1_fastp.json
    │   ├── Liver_ZT0_1_R1.fastq.gz
    │   ├── Liver_ZT0_1_R2.fastq.gz
    │   ├── Liver_ZT12_1_fastp_error.log
    │   ├── Liver_ZT12_1_fastp.html
    │   ├── Liver_ZT12_1_fastp.json
    │   ├── Liver_ZT12_1_R1.fastq.gz
    │   └── Liver_ZT12_1_R2.fastq.gz
    ├── Liver_ZT0_1_R1.fastq.gz
    ├── Liver_ZT0_1_R2.fastq.gz
    ├── Liver_ZT12_1_R1.fastq.gz
    ├── Liver_ZT12_1_R2.fastq.gz
    └── qcreports
        ├── Liver_ZT0_1_R1_fastqc.html
        ├── Liver_ZT0_1_R1_fastqc.zip
        ├── Liver_ZT0_1_R2_fastqc.html
        ├── Liver_ZT0_1_R2_fastqc.zip
        ├── Liver_ZT12_1_R1_fastqc.html
        ├── Liver_ZT12_1_R1_fastqc.zip
        ├── Liver_ZT12_1_R2_fastqc.html
        ├── Liver_ZT12_1_R2_fastqc.zip
        ├── multiqc_data
        │   ├── multiqc_citations.txt
        │   ├── multiqc_data.json
        │   ├── multiqc_fastqc.txt
        │   ├── multiqc_general_stats.txt
        │   ├── multiqc.log
        │   ├── multiqc_software_versions.txt
        │   └── multiqc_sources.txt
        └── multiqc_report.html


</details>

<details>
  <summary><i>Tasks</i></summary>
> Given: 
  
  <details>
    <summary><i>A. Task 1</i></summary>
> To run the pipeliine Note:
--- 
  </details>


  <details>
    <summary><i>B. Task 2</i></summary>

---
  </details>

</details>



***
