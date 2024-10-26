# Install the latest version of DEseq2
#Remember that you need to install the package only once but load it every time you start a new R session.
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)

##############################################################################
# Data prepare
# Load the counts filtered .txt file
data1 <- read.table("allfeaturecounts.tsv",row.names = 1, header = TRUE, sep = "\t")
dim(data1)
head(data1)

data <- data[, !(names(data) %in% c("Chr","Start" ,"End", "Strand","Length"))]
head(data)
dim(data)
# Subset counts file to protein_coding as gene_type if not done before then drop col gene_name/type
data <- data1 %>% filter(gene_type == "protein_coding")
head(data)
dim(data)
#,"gene_name", "gene_type"

data <- data[, !(names(data) %in% c("gene_type"))]


# Removing the singletons Note: it is better practice to filter at 0 now
# then filter your end results by their mean expression. But this is simpler*
data <- data[which(rowSums(data[, -1]) > 4), ]
dim(data)
colnames(data)<- c( "Genenames","Heart_ZT0_1","Heart_ZT0_2","Heart_ZT12_1","Heart_ZT12_2","Liver_ZT0_1","Liver_ZT0_2","Liver_ZT12_1","Liver_ZT12_2")

# Fix ensemble ids with version 1st by removing dot suffix 2nd by biomart and AnnotatinHub
rownames(data) <- gsub("\\..*", "", rownames(data))
write.csv(data, file = "proteincoding_geneids_name.csv")

#we can use thsi genenames to get gene name for id 
genenames <- data[,1,drop= FALSE]
genenames$ensemblids <- rownames(data)
genenames <- genenames[, c("ensemblids", names(genenames)[1])]
row.names(genenames)<- NULL
head(genenames)
write.csv(data, file = "geneids_name.csv")
#"gene_name", "gene_type"
data <- data[, !(names(data) %in% c("Genenames"))]
# RNA-seq data from two house mouse (Mus musculus) tissues (Heart, Liver) across 
# tow sampling times (ZT0, ZT12), with 1 biological replicates for each tissue and sampling time,
# resulting in a total of 16 paired-end FASTQ files.

org <- factor(c("m1","m2","m1","m2","m1","m2","m1","m2"))
typetissue <- factor(c("heart","heart","heart","heart","liver","liver","liver","liver"))
time <- factor(c("ZT0","ZT0","ZT12","ZT12","ZT0","ZT0","ZT12","ZT12"))
coldata <- data.frame(row.names = colnames(data), org, typetissue,time)
coldata
#coldata <- data.frame(row.names = colnames(data)[-1]
#################################################################################


ds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~typetissue + time + typetissue:time)
ds <- DESeq(ds)

vs <- vst(ds, blind = FALSE)
# Calculate the correlation matrix
cor_matrix <- cor(assay(ds))
# Visualize the correlation matrix as a heatmap
library(pheatmap)
pheatmap(cor_matrix, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         main = "Sample Correlation Heatmap")

pcaData <- plotPCA(vs,intgroup = c("typetissue","time"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

row_means <- rowMeans(normalized_counts)

ggplot(pcaData, aes(PC1, PC2, color=time, shape=typetissue)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#install.packages("GGally")
library(GGally)

#######
normalized_counts <- counts(ds, normalized = TRUE)

# Pairwise scatter plot matrix using base R
sample_data <- as.data.frame(normalized_counts[, c(1, 2, 5, 6)])

# Use ggpairs to create a scatter plot matrix
ggpairs(sample_data, 
        title = "4x4 Scatter for ZW0  Plot Matrix",
        lower = list(continuous = wrap("points", alpha = 0.5, size = 1)),
        upper = list(continuous = wrap("cor", size = 5))) # Use upper.panel = NULL to only show lower half



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
# Create a heatmap
pheatmap(heatmap_data, scale = "row", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", show_rownames = TRUE,
         main = "Top Variable Genes Heatmap")
head(heatmap_data)



# Dispersion is a measure of how much the expression of a gene varies across samples
plotDispEsts(ds)
