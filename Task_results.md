
<details open>
  <summary><i>Tasks</i></summary>
  
>**_Given:_**
> RNA-seq data from **two house mouse (Mus musculus) tissues (Heart, Liver)** across **two sampling times (ZT0, ZT12)**, with biological replicate for each tissue and sampling time, resulting in a total of **16 paired-end FASTQ files**.
> To analyze RNA-seq data, genome reference, and the corresponding GTF annotation file, with the analysis split into two main parts: __bioinformatics and statistical analysis__.
> Time taken for the pipeline to run without withouts 

  
  <details>
    <summary><i>A. Task 1</i></summary>
    
>Longest Substring Calculator
    
#### 1. a. Quality Control: 
- [X] [Fastqc](https://github.com/gunj007/RNA-Seq/tree/main/qcreports/rawfq_qc) Performed quality control using FastQC on rawfastq.gz 16
- [X] [MultiQC](https://github.com/gunj007/RNA-Seq/blob/main/qcreports/rawfq_qc/multiqc_report.html) Provided a summary report using MultiQC, adapter contamination seen. 
#### b. Adapter Trimming: 
- [X] using fastp. 
- [x] Provide a summary report detailing the percentage of reads trimmed. [Multiqc](https://github.com/gunj007/RNA-Seq/blob/main/qcreports/trimfq_qc/multiqc_report.html)
#### c. Genome Preparation: 
- [X] Prepared a [genome index using HISAT2](https://github.com/gunj007/RNA-Seq#genome-build--hisat2-genome-index) 
#### d. Alignment and Mapping: 
- [ ] Performed read alignment using HISAT2
- [ ] Provide alignment statistics and its visualization report, including the percentage of aligned reads, mapped reads, and potential issues with multi-mapping. 
#### e. Read Quantification: 
- [ ] Quantify gene expression using featureCounts or other suitable quantification software to generate a gene count expression matrix and provide associated statistical reports. 
- [ ] Output the results in a tabular format with protein coding genes ID as rows and samples as columns.



--- 
  </details>


  <details>
    <summary><i>B. Task 2</i></summary>
>Perform the differential expression and functional enrichment analysis using the gene count matrix created in bioinformatics analysis. 
    
#### a. Data reproducibility and pattern of variation: 
    
- [ ] Determine the reproducibility of the sample and then visualize the data using scatter plots and heatmaps. 
- [ ] Perform Principal Component Analysis (PCA) to investigate the overall pattern of variation across all samples in the dataset. 
#### b. Differential Expression Analysis: 
- [ ]  Use a factorial design in DESeq2 to model main effects (tissue and sampling time) and their interaction (formula: ~ tissue + time + tissue:time). Clearly report the statistical analysis of differentially expressed genes (DEGs) with 1) tissue-specific, 2) time-specific, and 3) interaction effects. 
- [ ]  Perform paired contrast analysis and detect DEGs between tissues at each sampling time. Visualize the results using volcano plots. 
- [ ]  Perform the clustering analysis for all DEGs and visualize the results using comprehensive heatmaps. 
- [ ]  For each statistical analysis group, display the expression patterns of the top DEGs and explain their biological significance. 
#### c. Functional Enrichment Analysis: 
- [ ]  For each statistical analysis group, perform GO/KEGG pathway enrichment analysis. Include bubble plots or bar plots for the top enriched pathways and provide tables summarizing key terms/pathways.

---
  </details>

</details>


>**_NOTE:_**  
> 1. To run the pipeliine on your system makesure you all the tools installed or refer [2.3 yml setup](https://github.com/gunj007/RNA-Seq?tab=readme-ov-file#23-installation-using-yml-rnaseq_envyml) and download the scripts/
> 2. In `counts.sh` change to `your_path_script/qc.sh` and for `hisat2.sh` before running
> 3. If genome is not built with the name genome then change it `your_genome_name` on line no. 

***
