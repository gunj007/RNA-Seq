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
  <summary><b>2. Installation Guide </b></summary>

  
 #### 2.1 System Info:
  - **System:** _Ubuntu 24.04 LTS_ `lsb_release -a`
  - **RAM - threads:** _7.45G - 11threads_  `htop`
  - **Specs:** 172G avail `df -h`

  ---
  
####  2.2 Conda Env Dependencies:

- Install [Miniconda](https://docs.anaconda.com/miniconda/)
- Add [Bioconda channels](https://bioconda.github.io/)
- Create env and INSTALLATION of TOOLS:
```batchfile
conda create -n ranaseq
conda activate rnaseq
```

```batchfile
# QC TOOLS
conda install bioconda::multiqc
conda install bioconda:fastqc
conda install -c bioconda fastp
```
 
```batchfile
# ALIGNMENT TOOLS
conda install bioconda::samtools
conda install bioconda::hisat2
conda install bioconda::subread
```

---
     
####  2.3 Installation using .yml: [rnaseq_env.yml](https://github.com/gunj007/RNA-Seq/tree/main/env)

```batchfile
conda env create -f rnaseq_env.yml
```

---

####  2.4 Tools versions: use tool_name `--version` or  `-v` to the version `--help` or `-h` for user guide of the tool.
     
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
  
####  3.1 Run commands:
  - FastqQC:
    - `-a ,  --adapters`    Specifies a non-default file which contains the list of  adapter sequences which will be explicity searched against the library. The file must contain sets of named adapters in the form name[tab]sequence.  Lines prefixed with a hash will be ignored.
    - `-q --quiet `      Suppress all progress messages on stdout and only report errors.
      
                    

```batchfile
fastqc -o output_dir *.gz --threads 10 -q
```

  - MultiQC:
    - **_ERR_**: fix try `conda update multiqc`  if error prevails config.py replace ymal.load to ymal.safe_load
      
            "/home/rgitbt/miniconda3/envs/multiqc/lib/python2.7/site-packages/multiqc-1.0.dev0-py2.7.egg/multiqc/utils/config.py:44:
           YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.

 

```batchfile
multiqc -o output_dir *zip
```

  - Fastp:
    -  `--adapter_sequence` the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto]
    -  `--adapter_sequence_r2` the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> eg ` --adapter_sequence=AGATCGGAAGAGC --adapter_sequence_r2=AGATCGGAAGAGC` if Multiqc report adapter is "Illumina Universal Adapter"
    - `--detect_adapter_for_pe` by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.
    - `--interleaved_in `  This option will result in interleaved FASTQ output for paired-end output. Disabled by default.


```batchfile
# "> /dev/null" as no --quiet option & 2>&1 as ran in to a ERR
fastp -i "$r1" -I "$r2" -o "$out_r1" -O "$out_r2" \
  -h "$html_report" -j "$json_report" > "$log_file" 2>&1 
```

  - Samtools:
```batchfile
samtools view -bS .sam > .bam
# you can also sort by coord incase of other aligners
```

  - Hisat2:
    - `--summary-file` <path> print alignment summary to this file.
    - `--time`
#### Genome Build : [Hisat2 Genome index](https://daehwankimlab.github.io/hisat2/howto/)
```batchfile
# Download Genome wget "link"
gunzip mgiGenome/GRCm39.primary_assembly.genome.fa.gz 
# BUID Genome this will create genome1.ht2 multiple files in genome directory
hisat2-build mgiGenome/GRCm39.primary_assembly.genome.fa genome
```

```batchfile
# RUN Command
# can use -p 10 for threads but requires more ram might crash ERR-137
hisat2 -x mgiGenome/genome/genome -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz -S trimmed.sam --quiet  --summary-file alignment_summary.txt --time
```

  - Subread(featureCounts):
```batchfile
featureCounts -p -t gene --extraAttributes gene_name,gene_type --primary -a annotation.gtf -o counts.txt 1.bam 2.bam nth.bam
```

  - Preprocessing: 
```batchfile
sed '1d' counts.txt > counts.tsv
# Subset on the basis of Protein_coding / exon etc 
```

---
  
####  3.2 Scripts:
  - QC Script: [qc.sh](https://github.com/gunj007/RNA-Seq/blob/main/scripts/qc.sh)
>This script performs _fastqc-multiqc-fastp-fastqc-multiqc_ for multiple files. Just provide input folder containing '.fast.gz', it will create a folder name `qcreports/` and it will create `qcreportstrim/` inside `fastp/` for trimmed reads

```batchfile
bash scripts/qc.sh ~/rawfastq
```

  - Alignment Script: [hisat2.sh](https://github.com/gunj007/RNA-Seq/blob/main/scripts/hisat2.sh)
>This scripts performs _hisat2-sam_to_bam_ provide trimmed fastq's folder ie. `fastp/` with the genome file and annotation file, it creates bam folder inside the input directory 
```batchfile
bash scripts/hisat2.sh ~/rawfastq ~/genome ~/annotations.gtf
```
---

####  3.3 Automation:
> To perform standard gene count matrix from raw FASTQ files run the [count.sh](https://github.com/gunj007/RNA-Seq/blob/main/scripts/count.sh)] script with the following command

```batchfile
#run "bash path_to_script_folder/count.sh path_to_rawfastq_folder/ path_to_genome_folder/ path_to_gtf-gff_file/.gtf
# main output of this script is featurecounts.tsv
bash scripts/count.sh ~/biostateai/raw_fastq ~/mgiGenome ~/mgiGenome/gencode.vM35.basic.annotation.gtf 
```
- Example Folder Structure:
```batchfile
../biostateai/
├── scripts <- Genome annotation file (.GTF/.GFF)
│   ├── count.sh
│   ├── hisat2.sh
│   └── qc.sh
├── raw_fastq
|   ├── bam
|   │   ├── all_bam.txt
|   │   ├── all_bam.txt.summary
|   |   ├── allfeaturecounts.tsv
|   │   ├── Liver_ZT0_1.bam
|   │   └── Liver_ZT12_1.bam
|   ├── fastp
|   │   ├── Liver_ZT0_1_fastp_error.log
|   │   ├── Liver_ZT0_1_fastp.html
|   │   ├── Liver_ZT0_1_fastp.json
|   │   ├── Liver_ZT0_1_R1.fastq.gz
|   │   ├── Liver_ZT0_1_R2.fastq.gz
|   │   ├── Liver_ZT12_1_fastp_error.log
|   │   ├── Liver_ZT12_1_fastp.html
|   │   ├── Liver_ZT12_1_fastp.json
|   │   ├── Liver_ZT12_1_R1.fastq.gz
|   │   └── Liver_ZT12_1_R2.fastq.gz
|   ├── Liver_ZT0_1_R1.fastq.gz
|   ├── Liver_ZT0_1_R2.fastq.gz
|   ├── Liver_ZT12_1_R1.fastq.gz
|   ├── Liver_ZT12_1_R2.fastq.gz
|   └── qcreports
|       ├── Liver_ZT0_1_R1_fastqc.html
|       ├── Liver_ZT0_1_R1_fastqc.zip
|       ├── Liver_ZT0_1_R2_fastqc.html
|       ├── Liver_ZT0_1_R2_fastqc.zip
|       ├── Liver_ZT12_1_R1_fastqc.html
|       ├── Liver_ZT12_1_R1_fastqc.zip
|       ├── Liver_ZT12_1_R2_fastqc.html
|       ├── Liver_ZT12_1_R2_fastqc.zip
|       ├── multiqc_data
|       │   ├── multiqc_citations.txt
|       │   ├── multiqc_data.json
|       │   ├── multiqc_fastqc.txt
|       │   ├── multiqc_general_stats.txt
|       │   ├── multiqc.log
|       │   ├── multiqc_software_versions.txt
|       │   └── multiqc_sources.txt
|       └── multiqc_report.html
└──mgiGenome/
    ├── gencode.vM35.basic.annotation.gtf
    ├── genome
    │   ├── genome.1.ht2
    │   ├── genome.2.ht2
    │   ├── genome.3.ht2
    │   ├── genome.4.ht2
    │   ├── genome.5.ht2
    │   ├── genome.6.ht2
    │   ├── genome.7.ht2
    │   └── genome.8.ht2
    └── GRCm39.primary_assembly.genome.fa

```

</details>

### [Task_results](https://github.com/gunj007/RNA-Seq/blob/main/Task_results)
---
