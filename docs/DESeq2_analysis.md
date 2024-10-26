# Differential Expression Analysis
>

#### load packages
```rscript
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)

```
#### pipeline output is required to be input here `allfeaturecounts.tsv`

```rscript
setwd("C:/Users/GUNJAN/Desktop/biostateAi/RNA-Seq/docs/pipeline_out")

> data1 <- read.table("allfeaturecounts.tsv",row.names = 1, header = TRUE, sep = "\t")
> dim(data1)
[1] 57132    15
> head(data1)
                      Chr   Start     End Strand
ENSMUSG00000102693.2 chr1 3143476 3144545      +
ENSMUSG00000064842.3 chr1 3172239 3172348      +
ENSMUSG00000051951.6 chr1 3276124 3741721      -
ENSMUSG00000102851.2 chr1 3322980 3323459      +
ENSMUSG00000103377.2 chr1 3435954 3438772      -
ENSMUSG00000104017.2 chr1 3445779 3448011      -
                     Length     gene_name
ENSMUSG00000102693.2   1070 4933401J01Rik
ENSMUSG00000064842.3    110       Gm26206
ENSMUSG00000051951.6 465598          Xkr4
ENSMUSG00000102851.2    480       Gm18956
ENSMUSG00000103377.2   2819       Gm37180
ENSMUSG00000104017.2   2233       Gm37363
                                gene_type
ENSMUSG00000102693.2                  TEC
ENSMUSG00000064842.3                snRNA
ENSMUSG00000051951.6       protein_coding
ENSMUSG00000102851.2 processed_pseudogene
ENSMUSG00000103377.2                  TEC
ENSMUSG00000104017.2                  TEC
                     rawfastq.bam.Heart_ZT0_1_sorted.bam
ENSMUSG00000102693.2                                   0
ENSMUSG00000064842.3                                   0
ENSMUSG00000051951.6                                  16
ENSMUSG00000102851.2                                   0
ENSMUSG00000103377.2                                   0
ENSMUSG00000104017.2                                   0
                     rawfastq.bam.Heart_ZT0_2_sorted.bam
ENSMUSG00000102693.2                                   0
ENSMUSG00000064842.3                                   0
ENSMUSG00000051951.6                                   9
ENSMUSG00000102851.2                                   0
ENSMUSG00000103377.2                                   0
ENSMUSG00000104017.2                                   0
                     rawfastq.bam.Heart_ZT12_1_sorted.bam
ENSMUSG00000102693.2                                    0
ENSMUSG00000064842.3                                    0
ENSMUSG00000051951.6                                   14
ENSMUSG00000102851.2                                    0
ENSMUSG00000103377.2                                    0
ENSMUSG00000104017.2                                    0
                     rawfastq.bam.Heart_ZT12_2_sorted.bam
ENSMUSG00000102693.2                                    0
ENSMUSG00000064842.3                                    0
ENSMUSG00000051951.6                                    3
ENSMUSG00000102851.2                                    0
ENSMUSG00000103377.2                                    0
ENSMUSG00000104017.2                                    0
                     rawfastq.bam.Liver_ZT0_1_sorted.bam
ENSMUSG00000102693.2                                   0
ENSMUSG00000064842.3                                   0
ENSMUSG00000051951.6                                   1
ENSMUSG00000102851.2                                   0
ENSMUSG00000103377.2                                   0
ENSMUSG00000104017.2                                   0
                     rawfastq.bam.Liver_ZT0_2_sorted.bam
ENSMUSG00000102693.2                                   0
ENSMUSG00000064842.3                                   0
ENSMUSG00000051951.6                                   0
ENSMUSG00000102851.2                                   0
ENSMUSG00000103377.2                                   0
ENSMUSG00000104017.2                                   0
                     rawfastq.bam.Liver_ZT12_1_sorted.bam
ENSMUSG00000102693.2                                    0
ENSMUSG00000064842.3                                    0
ENSMUSG00000051951.6                                    1
ENSMUSG00000102851.2                                    0
ENSMUSG00000103377.2                                    0
ENSMUSG00000104017.2                                    0
                     rawfastq.bam.Liver_ZT12_2_sorted.bam
ENSMUSG00000102693.2                                    0
ENSMUSG00000064842.3                                    0
ENSMUSG00000051951.6                                    0
ENSMUSG00000102851.2                                    0
ENSMUSG00000103377.2                                    0
ENSMUSG00000104017.2                                    0
```

###  tabular format with protein coding genes ID as rows and samples as columns.
> if added gene names can be dropped
- since in the script count file was not cleaned, 1st remove unnessecary columns, except gene_name and gene type,2nd rename colmns,3rd subset data to protein_coding as gene_type then drop gene_type column,
since there are many 0 counts will set sum to >4 as # Removing the singletons Note: it is better practice to filter at 0 now
# then filter your end results by their mean expression. But this is simpler*

```rscript
> data <- data1 %>% filter(gene_type == "protein_coding")
> dim(data)
[1] 21652    15

> data <- data[, !(names(data) %in% c("Chr","Start" ,"End", "Strand","Length"))]
> dim(data)
[1] 21652    10

> data <- data[, !(names(data) %in% c("gene_type"))]
> dim(data)
[1] 21652     9

> data <- data[which(rowSums(data[, -1]) > 4), ]
> dim(data)
[1] 15683     9

> colnames(data)
[1] "Genename"     "Heart_ZT0_1"  "Heart_ZT0_2"  "Heart_ZT12_1" "Heart_ZT12_2"
[6] "Liver_ZT0_1"  "Liver_ZT0_2"  "Liver_ZT12_1" "Liver_ZT12_2"

write.csv(data, file = "proteincoding_geneids_name.csv")
```
- since we have geneids with there version no.s we can either remove the version otherwise it will through a error or we can extract the gnene_name, gene id so while name we can label with gene name anyway is fine

```rscript
# rownames were 1st assigned to gene ids
> genenames <- data[,1,drop= FALSE]
> genenames$ensemblids <- rownames(data)
> genenames <- genenames[, c("ensemblids", names(genenames)[1])]
> row.names(genenames)<- NULL
> head(genenames)
             ensemblids Genename
1  ENSMUSG00000051951.6     Xkr4
2 ENSMUSG00000025900.14      Rp1
3 ENSMUSG00000025902.14    Sox17
4 ENSMUSG00000033845.14   Mrpl15
5  ENSMUSG00000104217.2  Gm37988
6 ENSMUSG00000033813.16    Tcea1

# drop gene_name
data <- data[, !(names(data) %in% c("Genename"))]
``

#### Deseq counts matrix

```rscript
> typetissue <- factor(c("heart","heart","heart","heart","liver","liver","liver","liver"))
> time <- factor(c("ZT0","ZT0","ZT12","ZT12","ZT0","ZT0","ZT12","ZT12"))
> coldata <- data.frame(row.names = colnames(data)[-1], org, typetissue,time)
> coldata
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

#### Visualize the correlation matrix as a heatmap, scatter & PCA
```rscript
# Apply variance-stabilizing transformation
vs <- vst(ds, blind = FALSE)
# Calculate the correlation matrix
cor_matrix <- cor(assay(ds))

######################2
library(pheatmap)
pheatmap(cor_matrix, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         main = "Sample Correlation Heatmap")



normalized_counts <- counts(ds, normalized = TRUE)
> colnames(normalized_counts)
[1] "Heart_ZT0_1"  "Heart_ZT0_2"  "Heart_ZT12_1" "Heart_ZT12_2" "Liver_ZT0_1" 
[6] "Liver_ZT0_2"  "Liver_ZT12_1" "Liver_ZT12_2"
#######################3
normalized_counts <- counts(ds, normalized = TRUE)
scatter_plot <- function(sample1, sample2, data) {
  plot(data[, sample1], data[, sample2],
       xlab = colnames(data)[sample1],
       ylab = colnames(data)[sample2],
       main = paste("Scatter Plot:", colnames(data)[sample1], "vs", colnames(data)[sample2]),
       pch = 20, col = rgb(0, 0, 1, 0.5))
  abline(lm(data[, sample2] ~ data[, sample1]), col = "red")
}

colnames(normalized_counts)
scatter_plot(4, 8, normalized_counts)



# Order genes by their means
normalized_counts_ordered <- normalized_counts[order(row_means, decreasing = TRUE), ]

################ Create a heatmap 4
pheatmap(normalized_counts_ordered,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",     # Scale rows to visualize variation
         show_rownames = FALSE)




pcaData <- plotPCA(vs,intgroup = c("typetissue","time"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time, shape=typetissue)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```
- 

```rscript






```


```rscript

```

```rscript

```
```rscript

```
```rscript

```
