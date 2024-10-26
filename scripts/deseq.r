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
############################################################################


# DESEQ2 object
ds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~typetissue + time + typetissue:time)
ds <- DESeq(ds)


# different comparisons available in DESeq
resultsNames(ds)

## will create function no repeat code
# Tissue-specific DEGs
restiss <- results(ds, name = "typetissue_liver_vs_heart")
head(restiss)
dim(restiss)

degtiss <- restiss[which(restiss$padj < 0.05), ]
summary(degtiss)
dim(degtiss)

tiss <- as.data.frame(degtiss)
head(tiss)
tiss$ensemblids <- rownames(tiss)
tiss <- tiss[,7]
head(tiss)
gene_names1 <- mapIds(org.Mm.eg.db,
                     keys = tiss,  # Use the cleaned column
                     column = "SYMBOL",    # Get gene symbols
                     keytype = "ENSEMBL",  # Key type is Ensembl
                     multiVals = "first")

head(gene_names1)
tiss <- as.data.frame(tiss)
type(gene_names1)

gene_names_only1 <- unname(gene_names1)

# Now add these gene names to your dataframe
tiss$GeneName <- gene_names_only1
sum(is.na(tiss$GeneName))

# Replace missing gene names with Ensembl IDs
tiss$GeneName[is.na(tiss$GeneName)] <- rownames(degtiss)[is.na(tiss$GeneName)]

# Now assign GeneName as rownames
rownames(degtiss) <- tiss$GeneName





ensemblIDs = rownames(restiss)
head(ensemblIDs)



symbols <- mapIds(org.Mm.eg.db, keys = ensemblIDs, keytype = "ENSEMBL", column="SYMBOL")
symbols = as.data.frame(symbols)
head(symbols)

# na gene we can extrac those gene na values from extracted gnen id and name col from bam file's 
na_rows <- symbols[is.na(symbols$symbols), ]

identical(rownames(symbols), rownames(restiss))

restiss = cbind.data.frame(symbols[,1], restiss)
colnames(restiss)[1] = "GeneSymbol"
dim(restiss)
write.csv(restiss, file = "DESeqResultstiss.csv")

#####################################################
# Time-specific DEGs
restime <- results(ds, name = "time_ZT12_vs_ZT0")
head(restime)
dim(restime)
degtime <- restime[which(restime$padj < 0.05), ]
summary(degtime)


tim <- as.data.frame(degtime)
head(tim)
tim$ensemblids <- rownames(tim)
tim <- tim[,7]
head(tim)
gene_names2 <- mapIds(org.Mm.eg.db,
                      keys = tim,  # Use the cleaned column
                      column = "SYMBOL",    # Get gene symbols
                      keytype = "ENSEMBL",  # Key type is Ensembl
                      multiVals = "first")

head(gene_names2)
tim <- as.data.frame(tim)
type(gene_names2)

gene_names_only2 <- unname(gene_names2)

# Now add these gene names to your dataframe
tim$GeneName <- gene_names_only2
sum(is.na(tim$GeneName))

# Replace missing gene names with Ensembl IDs
tim$GeneName[is.na(tim$GeneName)] <- rownames(degtime)[is.na(tim$GeneName)]

# Now assign GeneName as rownames
rownames(degtime) <- tim$GeneName




ensemblIDs = rownames(restime)
head(ensemblIDs)

symbols <- mapIds(org.Mm.eg.db, keys = ensemblIDs, keytype = "ENSEMBL", column="SYMBOL")
symbols = as.data.frame(symbols)
head(symbols)

# na gene we can extrac those gene na values from extracted gnen id and name col from bam file's 
na_rows <- symbols[is.na(symbols$symbols), ]

identical(rownames(symbols), rownames(restime))

restime = cbind.data.frame(symbols[,1], restime)
colnames(restime)[1] = "GeneSymbol"
dim(restime)
write.csv(restime, file = "DESeqResultstime.csv")

###################################################

# Interaction-specific DEGs (for example, heart at ZT0 vs liver at ZT0)
resint <- results(ds, name = "typetissueliver.timeZT12")
head(resint)
dim(resint)
# Filter DEGs with adjusted p-value < 0.05
degint <- resint[which(resint$padj < 0.05), ]
summary(degint)



int <- as.data.frame(degint)
head(int)
int$ensemblids <- rownames(int)
int <- int[,7]
head(int)
gene_names3 <- mapIds(org.Mm.eg.db,
                      keys = int,  # Use the cleaned column
                      column = "SYMBOL",    # Get gene symbols
                      keytype = "ENSEMBL",  # Key type is Ensembl
                      multiVals = "first")

head(gene_names3)
int <- as.data.frame(int)
type(gene_names3)

gene_names_only3 <- unname(gene_names3)

# Now add these gene names to your dataframe
int$GeneName <- gene_names_only3
sum(is.na(int$GeneName))

# Replace missing gene names with Ensembl IDs
int$GeneName[is.na(int$GeneName)] <- rownames(degint)[is.na(int$GeneName)]

# Now assign GeneName as rownames
rownames(degint) <- int$GeneName


ensemblIDs = rownames(resint)
head(ensemblIDs)

symbols <- mapIds(org.Mm.eg.db, keys = ensemblIDs, keytype = "ENSEMBL", column="SYMBOL")
symbols = as.data.frame(symbols)
head(symbols)

# na gene we can extrac those gene na values from extracted gnen id and name col from bam file's 
na_rows <- symbols[is.na(symbols$symbols), ]

identical(rownames(symbols), rownames(resint))

resint = cbind.data.frame(symbols[,1], resint)
colnames(resint)[1] = "GeneSymbol"
dim(resint)
write.csv(resint, file = "DESeqResultsint.csv")

#################################################################

# DEG analysis between tissues at ZT0

#BiocManager::install("EnhancedVolcano")
require(DOSE)
graphics.off()
library(EnhancedVolcano)
res0 <- results(ds, name = "typetissue_liver_vs_heart")
# Volcano plot for tiss--res/deg
# gene id ===rownames(res12) name ====res0$gene_symbol
EnhancedVolcano(res0,
                lab = res0$gene_symbol,
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


#gene name 
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = rownames(res0),  
                       column = "SYMBOL",      
                       keytype = "ENSEMBL",   
                       multiVals = "first")     

# Add gene symbols to your results data frame
res0$gene_symbol <- gene_symbols



# DEG analysis between tissues at ZT12
# resint--->degint
res12 <- results(ds, name ="typetissueliver.timeZT12")

EnhancedVolcano(res12,
                lab = (res12$gene_symbol),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'ZT12: Liver vs Heart',
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

#gene name 
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = rownames(res12),  
                       column = "SYMBOL",      
                       keytype = "ENSEMBL",   
                       multiVals = "first")     

# Add gene symbols to your results data frame
res12$gene_symbol <- gene_symbols

##############################################################

#do the same all groups tiss/0, int/12 , time0vs12

#res12 <- results(ds, name ="typetissueliver.timeZT12")
res012 <- results(ds, name = "time_ZT12_vs_ZT0")
# Filter significant DEGs tissue
degdata <- res012[which(res0$padj < 0.05), ]

# Extract normalized counts for these DEGs
normcounts <- assay(vst(ds))[rownames(degdata), ]

library(pheatmap)

# Generate a heatmap with clustering for rows (genes) and columns (samples)
pheatmap(normcounts,
         scale = "row",      
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("green", "white", "magenta"))(50),
         show_rownames = F,             
         main = "Clustering of All DEGs ZT12vs0")
# Top DEGs based on adjusted p-value
topdegs <- degdata[order(degdata$padj), ][1:10, ]
topgenes <- rownames(topdegs)

# Extract normalized counts for these top DEGs
topcounts <- normcounts[topgenes, ]
#genename
gene <- topdegs$gene_symbol
# Extract normalized counts for the top DEGs
topcounts <- normcounts[topgenes, ]

# Set row names to gene names
rownames(topcounts) <- gene
pheatmap(topcounts,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("green", "white", "magenta"))(50),
         main = "Top DEGs Expression Patterns z12vs0",
         show_rownames = TRUE)


#######################################
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
barplot(ego, showCategory = 10, title = "GO Biological Process Enrichment for Top DEGs zt12vs0")

############################################################################################