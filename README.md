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
|8. |RStudio|2022.12.0|
---
  </details>


<details open>
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
    - `--quiet`  print nothing to stderr except serious errors
    - `--met-file` send metrics to file at <path> (off)

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
    - `-T <int>`            Number of the threads. 1 by default.
    -  
    - ` -g <string>  `       Specify attribute type in GTF annotation. 'gene_id' by default. Meta-features used for read counting will be extracted from annotation using the provided value.
    -  ` --extraAttributes`   Extract extra attribute types from the provided GTF annotation and include them in the counting output. These attribute types will not be used to group features. If more than one attribute type is provided they should be separated by comma.


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
- Example Folder Structure: [link]()


</details>

<details open>
  <summary><b>4. RNAseq- Analysis RScripts </b></summary>

#### Tabular format genes & count matrix, reproducibility of the sample and then visualize the data using scatter plots and heatmaps

#### 4.1 Preprocessing the pipeline's counts.tsv
```Rscript
# Install the latest version of DEseq2
# the counts file from featurecouts cam be manipulated here to generate gene count matrix
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
# Load the counts filtered .txt file
data1 <- read.table("bam/allfeaturecounts.tsv",row.names = 1, header = TRUE, sep = "\t")
dim(data1)
[1] 57132    15

colnames(data1)
[1] "Chr"                                        
 [2] "Start"                                      
 [3] "End"                                        
 [4] "Strand"                                     
 [5] "Length"                                     
 [6] "gene_name"                                  
 [7] "gene_type"                                  
 [8] "...biostateai.rawfastq.bam.Heart_ZT0_1.bam" 
 [9] "...biostateai.rawfastq.bam.Heart_ZT0_2.bam" 
[10] "...biostateai.rawfastq.bam.Heart_ZT12_1.bam"
[11] "...biostateai.rawfastq.bam.Heart_ZT12_2.bam"
[12] "...biostateai.rawfastq.bam.Liver_ZT0_1.bam" 
[13] "...biostateai.rawfastq.bam.Liver_ZT0_2.bam" 
[14] "...biostateai.rawfastq.bam.Liver_ZT12_1.bam"
[15] "...biostateai.rawfastq.bam.Liver_ZT12_2.bam"

data <- data[, !(names(data) %in% c("Chr","Start" ,"End", "Strand","Length","gene_name", "gene_type"))]

#change colnames

# Subset counts file to protein_coding as gene_type if not done before then drop col gene_name/type
data <- data1 %>% filter(gene_type == "protein_coding")
dim(data)
[1] 21652    10

data <- data[which(rowSums(data[, -1]) > 4), ]
dim(data)
[1] 15496     8
#change colnames
colnames(data)<- c("Heart_ZT0_1","Heart_ZT0_2","Heart_ZT12_1","Heart_ZT12_2","Liver_ZT0_1","Liver_ZT0_2","Liver_ZT12_1","Liver_ZT12_2")

# Fix ensemble ids with version 1st by removing dot suffix 2nd by biomart and AnnotatinHub
rownames(data) <- gsub("\\..*", "", rownames(data))
write.csv(data, file = "rna_analysis/proteincoding_geneids.csv")

write.csv(data, file = "protein_coding_geneids.csv")
```

#### 4.2 geneids  to gene name
```
library(org.Mm.eg.db) 
# if required gene names can be achieved 
genenames <- data[,1,drop= FALSE]
genenames$ensemblids <- rownames(data)
genenames <- genenames[, c("ensemblids", names(genenames)[1])]
row.names(genenames)<- NULL
head(genenames)
```
- DESeq Object
```Rscript
# RNA-seq data from two house mouse (Mus musculus) tissues (Heart, Liver) across 
# tow sampling times (ZT0, ZT12), with 1 biological replicates for each tissue and sampling time,
# resulting in a total of 16 paired-end FASTQ files.

org <- factor(c("m1","m2","m1","m2","m1","m2","m1","m2"))
typetissue <- factor(c("heart","heart","heart","heart","liver","liver","liver","liver"))
time <- factor(c("ZT0","ZT0","ZT12","ZT12","ZT0","ZT0","ZT12","ZT12"))
coldata <- data.frame(row.names = colnames(data), org, typetissue,time)
coldata
             org typetissue time
Heart_ZT0_1   m1      heart  ZT0
Heart_ZT0_2   m2      heart  ZT0
Heart_ZT12_1  m1      heart ZT12
Heart_ZT12_2  m2      heart ZT12
Liver_ZT0_1   m1      liver  ZT0
Liver_ZT0_2   m2      liver  ZT0
Liver_ZT12_1  m1      liver ZT12
Liver_ZT12_2  m2      liver ZT12



> ds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~typetissue + time + typetissue:time)
> ds <- DESeq(ds)
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

``` 
#### Correlation
```Rscript
# Visualize the correlation matrix as a heatmap
library(pheatmap)
pheatmap(cor_matrix, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         main = "Sample Correlation Heatmap")
```
#### Top Variable genes
```Rscript
# Extract the top variable genes
library(matrixStats) 
top_genes <- head(order(rowVars(assay(ds), useNames = TRUE ), decreasing = TRUE), 20)
heatmap_data <- assay(ds)[top_genes, ]
#replace gene id to gene name
heat <- as.data.frame(heatmap_data)
head(heat)
heat$ensemblids <- rownames(heat)
heat <-heat[,9]
heat

library(org.Mm.eg.db) 
gene_names <- mapIds(org.Mm.eg.db,
                     keys = heat,  # Use the cleaned column
                     column = "SYMBOL",    # Get gene symbols
                     keytype = "ENSEMBL",  # Key type is Ensembl
                     multiVals = "first")

head(gene_names)
heat <- as.data.frame(heat)
type(gene_names)

gene_names_only <- unname(gene_names)

# Now add these gene names to your dataframe
heat$GeneName <- gene_names_only
rownames(heatmap_data) <- heat$GeneName
# Create a heatmap
pheatmap(heatmap_data, scale = "row", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", show_rownames = TRUE,
         main = "Top Variable Genes Heatmap")
head(heatmap_data)
```


#### PCA
```Rscript
pcaData <- plotPCA(vs,intgroup = c("typetissue","time"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time, shape=typetissue)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```

#### Differential Expression Analysis:

```Rscript

```
#### Volcano plot & Dot

```Rscript
library(EnhancedVolcano)

# Volcano plot for tiss--res/deg
EnhancedVolcano(degtiss,
                lab = rownames(degtiss),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Tissue-specific DEGs (Liver vs Heart)',
                pCutoff = 0.05,
                FCcutoff = 1.5,  # Set a fold change cutoff if needed
                pointSize = 2.0,
                labSize = 3.0)

# Create gene_type column
#degtiss$gene_type <- ifelse(degtiss$log2FoldChange > 0 & degtiss$padj < 0.05, "Upregulated",
#                            ifelse(degtiss$log2FoldChange < 0 & degtiss$padj < 0.05, "Downregulated", "Not DE"))

# res/degtime
EnhancedVolcano(degtime,
                lab = rownames(degtime),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Time-specific DEGs (ZT12 vs ZT0)',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 3.0)

# res/degint
EnhancedVolcano(resint,
                lab = rownames(resint),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Interaction-specific DEGs (Liver vs Heart at ZT12)',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 3.0)


require(DOSE)
# if plots are not previewing
graphics.off()

# change object gsetim
dotplot(gsetiss, showCategory=10, split=".sign") + facet_grid(.~.sign)+
  ggtitle("Top 10 GO Enrichment Terms (Tissue-Specific)")
)


```
#### GO
```Rscript

```
#### KEGG

```Rscript

```




  </details>

### [Task_results](https://github.com/gunj007/RNA-Seq/blob/main/Task_results.md)
---
