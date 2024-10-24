
<details>
  <summary><i>Tasks</i></summary>
  
>**_Given:_**
> RNA-seq data from **two house mouse (Mus musculus) tissues (Heart, Liver)** across **two sampling times (ZT0, ZT12)**, with biological replicate for each tissue and sampling time, resulting in a total of **16 paired-end FASTQ files**.
> To analyze RNA-seq data, genome reference, and the corresponding GTF annotation file, with the analysis split into two main parts: __bioinformatics and statistical analysis__.
> Time taken for the pipeline to run without withouts 

  
  <details>
    <summary><i>A. Task 1</i></summary>
    
>Longest Substring Calculator
    
#### 1. a. Quality Control: 
- [X] [Fastqc](https://github.com/gunj007/RNA-Seq/tree/main/qcreports/fastqc_out) Perform quality control using FastQC 
- [X] [MultiQC](https://github.com/gunj007/RNA-Seq/blob/main/qcreports/multiqc_out/multiqc_report.html) Provide a summary report using MultiQC with key quality metrics such as sequence count per sample, per base sequence quality, read length distribution, sequence duplication level and any adapter contamination. 
#### b. Adapter Trimming: 
- [X] fastp or other suitable preprocessing and quality control software. 
- [x] Provide a summary report detailing the percentage of reads trimmed and any improvement in quality metrics after trimming. 
#### c. Genome Preparation: 
- [X] Prepare a genome index using HISAT2 or other suitable alignment software for the provided genome reference. 
#### d. Alignment and Mapping: 
- [X] Perform read alignment using HISAT2 or other suitable alignment software. 
- [X] Provide alignment statistics and its visualization report, including the percentage of aligned reads, mapped reads, and potential issues with multi-mapping. 
#### e. Read Quantification: 
i. Quantify gene expression using featureCounts or other suitable quantification software to generate a gene count expression matrix and provide associated statistical reports. 
ii. Output the results in a tabular format with protein coding genes ID as rows and samples as columns.





>**_NOTE:_**  
> 1. To run the pipeliine on your system makesure you all the tools installed or refer 2.3 and download the scripts/
> 2. In `counts.sh` change to `your_path_script/qc.sh` and for `hisat2.sh` before running
> 3. If genome is not built with the name genome then change it `your_genome_name` on line no. 

--- 
  </details>


  <details>
    <summary><i>B. Task 2</i></summary>

---
  </details>

</details>



***
