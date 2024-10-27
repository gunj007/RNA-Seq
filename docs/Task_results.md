

## RAN-Seq END TO END Analysis Index

<details>
  <summary><i>INFO</i></summary>
  
>**_Given:_**
> RNA-seq data from **two house mouse (Mus musculus) tissues (Heart, Liver)** across **two sampling times (ZT0, ZT12)**, with biological replicate for each tissue and sampling time, resulting in a total of **16 paired-end FASTQ files**.
> To analyze RNA-seq data, genome reference, and the corresponding GTF annotation file, with the analysis split into two main parts: __bioinformatics and statistical analysis__.
> Time taken for the pipeline [count.sh](https://github.com/gunj007/RNA-Seq/blob/main/scripts/count.sh) to run ~5hrs with ~7G RAM 127G Storage

  </details>
  
  <details open>
    <summary><i>Task Check list</i></summary>

> __**NOTE**__
> Task 2 [Report.pdf](https://github.com/gunj007/RNA-Seq/blob/main/docs/Report.pdf)
> Copy .html this link and open in [htmlpreview](https://htmlpreview.github.io/?) top view the reports


| Sr. Check | INDEX | TOOLS/CODE/EXPLANATION/SCRIPT LINKS | RESULTS |
|:---------|:------:|:---------------------------:|:---------:|
|<ul><li>- [x] 1.</li> | Quality Control | [FASTQC+MULTIQC](https://github.com/gunj007/RNA-Seq/blob/main/README.md#311-qc), [qc.sh](https://github.com/gunj007/RNA-Seq/blob/main/scripts/qc.sh) | Before Trim: [FastqcReports](https://github.com/gunj007/RNA-Seq/tree/main/qcreports/rawfq_qc), [MultiQC](https://github.com/gunj007/RNA-Seq/blob/main/docs/pipeline_out/qcreports/multiqc_report.html) -> adapter contamination seen in Reverse reads sequence(R2); [POST TRIM MULTIQC](https://github.com/gunj007/RNA-Seq/blob/main/docs/pipeline_out/qcreports/fastprepo/multiqc_report.html)|
|<ul><li>- [x] 2.</li> | Adapter Trimming | [FASTP](https://github.com/gunj007/RNA-Seq/blob/main/README.md#312-trimming) | [Multiqc](https://github.com/gunj007/RNA-Seq/blob/main/docs/pipeline_out/qcreports/fastprepo/multiqc_report.html), [fastp summary](https://github.com/gunj007/RNA-Seq/tree/main/docs/pipeline_out/fastprepo)
|<ul><li>- [x] 3.</li> | Genome Preparation | [HISAT2](https://github.com/gunj007/RNA-Seq/blob/main/README.md#313-genome-build--hisat2-genome-index) | |
|<ul><li>- [x] 4.</li> | Alignment and Mapping | [HISAT2](https://github.com/gunj007/RNA-Seq/blob/main/README.md#314-alignment--mapping), [hisat2.sh](https://github.com/gunj007/RNA-Seq/blob/main/scripts/hisat2.sh), [flagsat.r](https://github.com/gunj007/RNA-Seq/blob/main/scripts/flagstat.r) | [Read alignment reports](https://github.com/gunj007/RNA-Seq/tree/main/bin/samsum6) ,[Alignment statistics](https://github.com/gunj007/RNA-Seq/blob/main/docs/supplementary/1flagstatAlignmentSummary.csv) and its [visualization report](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/1Rplotflagstat_alignment.jpeg) |
|<ul><li>- [x] 5.</li> | Read Quantification | [featureCounts](https://github.com/gunj007/RNA-Seq/blob/main/README.md#315-read-quantification) ,[count.sh](https://github.com/gunj007/RNA-Seq/blob/main/scripts/count.sh) | [Protein_coding_genes.csv](https://github.com/gunj007/RNA-Seq/blob/main/docs/supplementary/proteincoding_geneids_name.csv) |
|<ul><li>- [x] 6.</li> | Data reproducibility and pattern of variation | [Data variation plot](https://github.com/gunj007/RNA-Seq/blob/main/docs/DESeq2_analysis.md#data-variation) [datavariation.r](https://github.com/gunj007/RNA-Seq/blob/main/scripts/datavariation.r) |  [Scatter plot](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/3Rplotscatterhvsl12.jpeg) & [heatmaps- orrelation matrix](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/2Rplotcor_matrix.jpeg), [Top variable genes](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/6RplotTopvariablegenesheat.jpeg), [Principal Component Analysis (PCA)](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/5RplotPCAnormal.jpeg) |
|<ul><li>- [x] 7.</li> | Differential Expression Analysis | [DESEQ2](https://github.com/gunj007/RNA-Seq/blob/main/docs/DESeq2_analysis.md#differential-expression-analysis-1); [DEG tissue vs time](https://github.com/gunj007/RNA-Seq/blob/main/docs/DESeq2_analysis.md#perform-paired-contrast-analysis-and-detect-degs-between-tissues-at-each-sampling-time-visualize-the-results-using-volcano-plots); [Expression Patterns and there biological significance](https://github.com/gunj007/RNA-Seq/blob/main/docs/DESeq2_analysis.md#heatmaps-top10degs-display-the-expression-patterns-of-the-top-degs-and-explain-their-biological-significance) [deseq.r](https://github.com/gunj007/RNA-Seq/blob/main/scripts/deseq.r) | DEG: 1) tissue-specific : [genes.csv](https://github.com/gunj007/RNA-Seq/blob/main/docs/supplementary/DESeqResultstiss.csv), 2) time-specific : [genes.csv](https://github.com/gunj007/RNA-Seq/blob/main/docs/supplementary/DESeqResultstime.csv), 3) interaction effects : [gene.csv](https://github.com/gunj007/RNA-Seq/blob/main/docs/supplementary/DESeqResultsint.csv); Time vs Tissue: [volcanoplotZT0](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/7Rplotdispersion.jpeg) AND [volcanoplotZT12](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/9Rplotvolcanonam12.jpeg) ; Top 10 DEGs [1st group](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/10tiss0Rplottopdegnam.jpeg) ,[2nd group](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/10tint1212Rplottop10degnamheat.jpeg) , [3rd group](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/10timeRplottop10genesnam.jpeg): Top 10 biological significance bar plots: [1st group](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/11tiss0Rplotbargotissz0.jpeg), [2nd](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/11timeRplotbar.jpeg), [3rd](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/11int12Rplotbargo.jpeg)  |
|<ul><li>- [x] 8.</li> | Functional Enrichment Analysis | [ENRICHGO](https://github.com/gunj007/RNA-Seq/blob/main/docs/Function_enrichment_analysis.md), [functional](https://github.com/gunj007/RNA-Seq/blob/main/scripts/functional.r)| Top enriched bar plot: enrich GO summary - [1st](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/12tissRplot10go.png): [summary](https://github.com/gunj007/RNA-Seq/blob/main/docs/supplementary/go_enrichment_tiss.csv), [2nd](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/12timeRplot10go.png):  [summary](https://github.com/gunj007/RNA-Seq/blob/main/docs/supplementary/go_enrichment_time.csv), [3rd](https://github.com/gunj007/RNA-Seq/blob/main/docs/plots/12intRplot10go.png) : [summary](https://github.com/gunj007/RNA-Seq/blob/main/docs/supplementary/go_enrichment_int.csv) |


---
  

</details>


>**_NOTE:_**  
> 1. To run the pipeliine on your system makesure you all the tools installed or refer [2. Installation Guide](https://github.com/gunj007/RNA-Seq/blob/main/README.md) and download the scripts/
> 2. In `counts.sh` change to `your_path_script/qc.sh` and for `hisat2.sh` before running
> 3. If genome is not built with the name genome then change it `your_genome_name` on line no. 

***
