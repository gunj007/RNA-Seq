## RNA-Seq Pipeline and Analysis

<details open>
  <summary><b>Introduction</b></summary>
  
  RNA-Seq is a sequencing method used to determine gene expression levels. [RNA-Seq](https://pmc.ncbi.nlm.nih.gov/articles/PMC6096346/) data originates from extracted RNA that was reverse transcribed into DNA.
There are many readyly available RNA-seq pipeline workflow eg. [nf-core](https://nf-co.re/rnaseq/3.14.0/), [PiGx](https://bioinformatics.mdc-berlin.de/pigx_docs/pigx-rna-seq.html), etc. [Galaxy](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html) is a popular framework software for constructing and managing the execution of workflows. 
  
</details>

<details>
  <summary><b>1.1 Tools & Workflow </b></summary>

1. Quality Control (QC): [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://github.com/MultiQC/MultiQC), & [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
2. Improving Quality (Trimming and Filtering): [fastp](https://github.com/OpenGene/fastp), [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic), [Cutadapt](https://cutadapt.readthedocs.io/en/stable/), & [Skewer](https://github.com/relipmoc/skewer)
3. Read Alignment: [HISAT2](https://daehwankimlab.github.io/hisat2/manual/), [STAR](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html), [TopHat2](https://ccb.jhu.edu/software/tophat/manual.shtml), [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
4. Transcript Assembly and Quantification: [StringTie](), [Cufflinks](), [Salmon](), [Kallisto]()
5. Differential Expression Analysis: [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), enrichplot, pathview, edgeR, limma- R packages
6. Gene Ontology (GO) and Pathway Analysis: [GSEA (Gene Set Enrichment Analysis)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page), clusterProfiler, EnhancedVolcano, DAVID
7. Visualization: [IGV (Integrative Genomics Viewer)](https://igv.org/doc/desktop/), pheatmap, ggplot2(R packages)
8. Single-Cell RNA-seq Specific Tools: CellRanger, Seurat, & SC3

</details>
  
<details>
  <summary><b>1.2 Installation Guide</b></summary>
  
  1. Conda Env:

     
  3. yml: Download 

```
conda env create -f rnaseq_env.yml
```

  4. System Info:


  5. Tools --versions:
     
|Sr.no|Tools|Version|
|:----|:----|:-----:|
|1. |Miniconda|conda 24.5.0|
|2. |Python|Python 3.12.4|
|3. |FastQC|FastQC v0.12.1|
|4. |MultiQC|version 1.18|
|5. |HISAT2|version 2.2.1|
|6. |samtools|samtools 1.19.2|
|7. |Subread|featureCounts v2.0.1|


  </details>


<details>
  <summary><b>1.3 Pipeline </b></summary>
  1 . Run cmds: give flg links
  2. R script


  <details>
    <summary><i>A. Task 1</i></summary>
  
  </details>


  <details>
    <summary><i>B. Task 2</i></summary>
  
  </details>

</details>
