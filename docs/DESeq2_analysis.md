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
> head(data1)                      Chr   Start     End Strand Length     gene_name            gene_type rawfastq.bam.Heart_ZT0_1_sorted.bam
ENSMUSG00000102693.2 chr1 3143476 3144545      +   1070 4933401J01Rik                  TEC                                   0
ENSMUSG00000064842.3 chr1 3172239 3172348      +    110       Gm26206                snRNA                                   0
ENSMUSG00000051951.6 chr1 3276124 3741721      - 465598          Xkr4       protein_coding                                  16
ENSMUSG00000102851.2 chr1 3322980 3323459      +    480       Gm18956 processed_pseudogene                                   0
ENSMUSG00000103377.2 chr1 3435954 3438772      -   2819       Gm37180                  TEC                                   0
ENSMUSG00000104017.2 chr1 3445779 3448011      -   2233       Gm37363                  TEC                                   0
                     rawfastq.bam.Heart_ZT0_2_sorted.bam rawfastq.bam.Heart_ZT12_1_sorted.bam rawfastq.bam.Heart_ZT12_2_sorted.bam
ENSMUSG00000102693.2                                   0                                    0                                    0
ENSMUSG00000064842.3                                   0                                    0                                    0
ENSMUSG00000051951.6                                   9                                   14                                    3
ENSMUSG00000102851.2                                   0                                    0                                    0
ENSMUSG00000103377.2                                   0                                    0                                    0
ENSMUSG00000104017.2                                   0                                    0                                    0
                     rawfastq.bam.Liver_ZT0_1_sorted.bam rawfastq.bam.Liver_ZT0_2_sorted.bam rawfastq.bam.Liver_ZT12_1_sorted.bam
ENSMUSG00000102693.2                                   0                                   0                                    0
ENSMUSG00000064842.3                                   0                                   0                                    0
ENSMUSG00000051951.6                                   1                                   0                                    1
ENSMUSG00000102851.2                                   0                                   0                                    0
ENSMUSG00000103377.2                                   0                                   0                                    0
ENSMUSG00000104017.2                                   0                                   0                                    0
                     rawfastq.bam.Liver_ZT12_2_sorted.bam
ENSMUSG00000102693.2                                    0
ENSMUSG00000064842.3                                    0
ENSMUSG00000051951.6                                    0
ENSMUSG00000102851.2                                    0
ENSMUSG00000103377.2                                    0
ENSMUSG00000104017.2                                    0
> 
```

###  tabular format with protein coding genes ID as rows and samples as columns.
> if added gene names can be dropped
- since in the script count file was not cleaned, 1st remove unnessecary columns, except gene_name and gene type,2nd rename colmns,3rd subset data to protein_coding as gene_type then drop gene_type column,
since there are many 0 counts will set sum to >4 as # Removing the singletons
> Note: it is better practice to filter at 0 now
then filter your end results by their mean expression. 

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
#######################3plot
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
normalized_counts_ordered <- 


pcaData <- plotPCA(vs,intgroup = c("typetissue","time"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
########333 5vsPCA
ggplot(pcaData, aes(PC1, PC2, color=time, shape=typetissue)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```
-  zvbfgfgb

```rscript

normalized_counts[order(row_means, decreasing = TRUE), ]

################ Create a heatmap 4
# Order genes by their means
normalized_counts_ordered <- normalized_counts[order(row_means, decreasing = TRUE), ]

# Create a heatmap genes
pheatmap(normalized_counts_ordered,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",     # Scale rows to visualize variation
         show_rownames = FALSE)
# Perform PCA on the normalized counts
pca_result <- prcomp(t(normalized_counts), center = TRUE, scale. = TRUE)
pca_variance <- (pca_result$sdev)^2
pca_variance_explained <- pca_variance / sum(pca_variance) 
pca_data <- as.data.frame(pca_result$x)
pca_data$condition <- colData(ds)$org 

ggplot(pca_data, aes(PC1, PC2, color = time, shape = typetissue)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() +
  ggtitle("PCA of Samples by Tissue Type and Time") +
  theme_minimal() +  # Optional: A cleaner theme
  theme(legend.position = "top")
#####
```

- dgfbsxfvsdv top variablr genes gene id to genename


```rscript
# Extract the top variable genes
library(matrixStats) 

library(org.Mm.eg.db) 
top_genes <- head(order(rowVars(assay(ds), useNames = TRUE ), decreasing = TRUE), 20)
heatmap_data <- assay(ds)[top_genes, ]
#replace gene id to gene name
heat <- as.data.frame(heatmap_data)
head(heat)
heat$ensemblids <- rownames(heat)
heat <-heat[,9]
heat

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
# Create a heatmap ############6
pheatmap(heatmap_data, scale = "row", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", show_rownames = TRUE,
         main = "Top Variable Genes Heatmap")
head(heatmap_data)

###################7
# Dispersion is a measure of how much the expression of a gene varies across samples
plotDispEsts(ds)


```

#### Differential Expression Analysis: 
i. Use a factorial design in DESeq2 to model main effects (tissue and sampling time) and their interaction (formula: ~ tissue + time + tissue:time). Clearly report the statistical analysis of differentially expressed genes (DEGs) with 1) tissue-specific, 2) time-specific, and 3) interaction effects. 


```rscript
> resultsNames(ds)
[1] "Intercept"                 "typetissue_liver_vs_heart"
[3] "time_ZT12_vs_ZT0"          "typetissueliver.timeZT12" 
```

- now we get 3names from results below deg for tissue change the name and o the same other 2 "time_ZT12_vs_ZT0"          "typetissueliver.timeZT12"

```rscript
> # Tissue-specific DEGs
> restiss <- results(ds, name = "typetissue_liver_vs_heart")

> dim(restiss)
[1] 15683     6

```


```rscript
> degtiss <- restiss[which(restiss$padj < 0.05), ]
> summary(degtiss)

out of 9822 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4851, 49%
LFC < 0 (down)     : 4971, 51%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> dim(degtiss)
[1] 9822    6
> 
> tiss <- as.data.frame(degtiss)
> head(tiss)
                      baseMean log2FoldChange     lfcSE      stat
ENSMUSG00000025900   12.924588     -7.1064968 1.6380102 -4.338494
ENSMUSG00000025902   88.150219     -3.0825541 0.3860948 -7.983931
ENSMUSG00000033845 1082.326791      0.4769393 0.1568420  3.040890
ENSMUSG00000104217    7.292742      2.6575155 1.0915134  2.434707
ENSMUSG00000033813  490.642842      1.1889130 0.2087431  5.695581
ENSMUSG00000033793  813.568758      1.5389819 0.1685540  9.130500
                         pvalue         padj
ENSMUSG00000025900 1.434626e-05 3.313754e-05
ENSMUSG00000025902 1.417449e-15 6.200414e-15
ENSMUSG00000033845 2.358802e-03 4.257253e-03
ENSMUSG00000104217 1.490385e-02 2.421955e-02
ENSMUSG00000033813 1.229525e-08 3.619896e-08
ENSMUSG00000033793 6.818369e-20 3.571233e-19
> tiss$ensemblids <- rownames(tiss)
> tiss <- tiss[,7]
> head(tiss)
[1] "ENSMUSG00000025900" "ENSMUSG00000025902" "ENSMUSG00000033845"
[4] "ENSMUSG00000104217" "ENSMUSG00000033813" "ENSMUSG00000033793"
> gene_names1 <- mapIds(org.Mm.eg.db,
+                      keys = tiss,  # Use the cleaned column
+                      column = "SYMBOL",    # Get gene symbols
+                      keytype = "ENSEMBL",  # Key type is Ensembl
+                      multiVals = "first")
'select()' returned 1:many mapping between keys and columns
> 
> head(gene_names1)
ENSMUSG00000025900 ENSMUSG00000025902 ENSMUSG00000033845 
             "Rp1"            "Sox17"           "Mrpl15" 
ENSMUSG00000104217 ENSMUSG00000033813 ENSMUSG00000033793 
                NA            "Tcea1"          "Atp6v1h" 
> tiss <- as.data.frame(tiss)
> type(gene_names1)
[1] "character"
> 
> gene_names_only1 <- unname(gene_names1)
> 
> # Now add these gene names to your dataframe
> tiss$GeneName <- gene_names_only1
> sum(is.na(tiss$GeneName))
[1] 32
> 
> # Replace missing gene names with Ensembl IDs
> tiss$GeneName[is.na(tiss$GeneName)] <- rownames(degtiss)[is.na(tiss$GeneName)]
> 
> # Now assign GeneName as rownames
> rownames(degtiss) <- tiss$GeneName
> 
> 
> 
> 
> 
> ensemblIDs = rownames(restiss)
> head(ensemblIDs)
[1] "ENSMUSG00000051951" "ENSMUSG00000025900" "ENSMUSG00000025902"
[4] "ENSMUSG00000033845" "ENSMUSG00000104217" "ENSMUSG00000033813"
> 
> 
> 
> symbols <- mapIds(org.Mm.eg.db, keys = ensemblIDs, keytype = "ENSEMBL", column="SYMBOL")
'select()' returned 1:many mapping between keys and columns
> symbols = as.data.frame(symbols)
> head(symbols)
                   symbols
ENSMUSG00000051951    Xkr4
ENSMUSG00000025900     Rp1
ENSMUSG00000025902   Sox17
ENSMUSG00000033845  Mrpl15
ENSMUSG00000104217    <NA>
ENSMUSG00000033813   Tcea1
> 
> # na gene we can extrac those gene na values from extracted gnen id and name col from bam file's 
> na_rows <- symbols[is.na(symbols$symbols), ]
> 
> identical(rownames(symbols), rownames(restiss))
[1] TRUE
> 
> restiss = cbind.data.frame(symbols[,1], restiss)
> colnames(restiss)[1] = "GeneSymbol"
> dim(restiss)
[1] 15683     7
> write.csv(restiss, file = "DESeqResultstiss.csv")
```

```rscript

dim(restiss)
[1] 15683     7
Warning messages:
1: In doTryCatch(return(expr), name, parentenv, handler) :
  display list redraw incomplete
2: In doTryCatch(return(expr), name, parentenv, handler) :
  invalid graphics state
3: In doTryCatch(return(expr), name, parentenv, handler) :
  invalid graphics state
> dim(restime)
[1] 15683     7
> dim(resint)
[1] 15683     7
> head(restiss)
                   GeneSymbol    baseMean log2FoldChange     lfcSE      stat       pvalue         padj
ENSMUSG00000051951       Xkr4    3.610031     -3.8078049 1.8766466 -2.029047 4.245346e-02 6.372669e-02
ENSMUSG00000025900        Rp1   12.924588     -7.1064968 1.6380102 -4.338494 1.434626e-05 3.313754e-05
ENSMUSG00000025902      Sox17   88.150219     -3.0825541 0.3860948 -7.983931 1.417449e-15 6.200414e-15
ENSMUSG00000033845     Mrpl15 1082.326791      0.4769393 0.1568420  3.040890 2.358802e-03 4.257253e-03
ENSMUSG00000104217       <NA>    7.292742      2.6575155 1.0915134  2.434707 1.490385e-02 2.421955e-02
ENSMUSG00000033813      Tcea1  490.642842      1.1889130 0.2087431  5.695581 1.229525e-08 3.619896e-08
> 
> head(resint)
                   GeneSymbol    baseMean log2FoldChange     lfcSE        stat     pvalue      padj
ENSMUSG00000051951       Xkr4    3.610031     1.12647760 2.6840980  0.41968571 0.67471506        NA
ENSMUSG00000025900        Rp1   12.924588     0.89834163 2.3196093  0.38728144 0.69854787        NA
ENSMUSG00000025902      Sox17   88.150219    -1.48053099 0.7049154 -2.10029607 0.03570280 0.1932835
ENSMUSG00000033845     Mrpl15 1082.326791     0.48534227 0.2243601  2.16322909 0.03052356 0.1750813
ENSMUSG00000104217       <NA>    7.292742     0.03868936 1.6627975  0.02326763 0.98143679        NA
ENSMUSG00000033813      Tcea1  490.642842     0.20525830 0.2992159  0.68598732 0.49272111 0.7695700
> 
> head(restime)
                   GeneSymbol    baseMean log2FoldChange     lfcSE        stat    pvalue      padj
ENSMUSG00000051951       Xkr4    3.610031    -0.24148640 1.3518338 -0.17863616 0.8582234        NA
ENSMUSG00000025900        Rp1   12.924588     0.05591804 0.8002888  0.06987232 0.9442953        NA
ENSMUSG00000025902      Sox17   88.150219     0.32307761 0.3078181  1.04957318 0.2939144 0.6962510
ENSMUSG00000033845     Mrpl15 1082.326791    -0.20002047 0.1583437 -1.26320435 0.2065158 0.6041844
ENSMUSG00000104217       <NA>    7.292742    -0.44490342 1.2862438 -0.34589354 0.7294227        NA
ENSMUSG00000033813      Tcea1  490.642842     0.02618329 0.2131760  0.12282477 0.9022459 0.9763164
> 
```







#### Perform paired contrast analysis and detect DEGs between tissues at each sampling time. Visualize the results using volcano plots. 

```rscript
library(EnhancedVolcano)
res0 <- results(ds, name = "typetissue_liver_vs_heart")
# Volcano plot for tiss--res/deg 
# DEG analysis between tissues at ZT0 ZT12
# gene id ===rownames(res12) name ====res0$gene_symbol
EnhancedVolcano(res0,
                lab = rownames(res0),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'ZT0: Liver vs Heart',
                pCutoff = 0.05,
                FCcutoff = 1,  # Set log2 fold change threshold for highlighting
                pointSize = 2,
                labSize = 3,
                col = c("grey30", "forestgreen", "royalblue", "red2"), # Custom colors
                colAlpha = 0.8,  # Transparency
                legendLabels = c("NS", "Log2FC > 1", "p < 0.05", "Upregulated or Downregulated"),
                legendPosition = "top",
                xlab = "Log2 Fold Change",
                ylab = "-Log10 Adjusted P-value")
```
```rscript
#gene name 
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = rownames(res0),  
                       column = "SYMBOL",      
                       keytype = "ENSEMBL",   
                       multiVals = "first")     

# Add gene symbols to your results data frame
res0$gene_symbol <- gene_symbols

```
#### heatmaps top10degs display the expression patterns of the top DEGs and explain their biological significance. 
do the same for all 3 groups 
```rscript
#do the same all groups tiss/0, int/12 , time0vs12

#res12 <- results(ds, name ="typetissueliver.timeZT12")
res012 <- results(ds, name = "time_ZT12_vs_ZT0")
# Filter significant DEGs tissue
degdata <- res0[which(res0$padj < 0.05), ]

# Extract normalized counts for these DEGs
normcounts <- assay(vst(ds))[rownames(degdata), ]

library(pheatmap)

# Generate a heatmap with clustering for rows (genes) and columns (samples)
pheatmap(normcounts,
         scale = "row",      
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("navy", "white", "firebrick"))(50),
         show_rownames = FALSE,             
         main = "Clustering of All DEGs")
# Top DEGs based on adjusted p-value
topdegs <- degdata[order(degdata$padj), ][1:10, ]
topgenes <- rownames(topdegs)

# Extract normalized counts for these top DEGs
topcounts <- normcounts[topgenes, ]

pheatmap(topcounts,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("navy", "white", "firebrick"))(50),
         main = "Top DEGs Expression Patterns")


library(clusterProfiler)
library(org.Mm.eg.db)

# Run GO enrichment for the top DEGs
ego <- enrichGO(gene = topgenes,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",                  # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

# View the most significant biological processes
barplot(ego, showCategory = 10, title = "GO Biological Process Enrichment for Top DEGs")
```



***