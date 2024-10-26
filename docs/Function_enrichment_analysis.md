# Function Enrichment Analysis
> Enrichment analysis in functional genomics helps identify over-represented biological functions or pathways within a list of genes differentially expressed genes, interpreting biological significance.
[](https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html)

```rscript
#install and load the required packages
#BiocManager::install("clusterProfiler"")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")

library(enrichplot)
library(ggplot2)
library(clusterProfiler)
library(pathview)
library(org.Mm.eg.db)

res0 <- results(ds, name ="typetissue_liver_vs_heart")
res12 <- results(ds, name ="typetissueliver.timeZT12")
res012 <- results(ds, name = "time_ZT12_vs_ZT0")

require(DOSE)
graphics.off()
```

#### 1. Gene Ontology (GO)

```rscript
# you read csv or load object
res = read.csv("DESeqResults.csv")
deg_genes <- rownames(res12[!is.na(res12$padj) & res12$padj < 0.05, ])
degenes <- mapIds(org.Mm.eg.db,
                            keys = deg_genes,
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")

# Remove any NA values from unmapped genes
degenes <- na.omit(degenes)

# GO Enrichment Analysis
go_results <- enrichGO(gene = degenes,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL", 
                       ont = "BP",          # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# Bubble plot for GO enrichment results
png(filename = "12intRplot10go.png", width = 800, height = 700, res = 100)
dotplot(go_results, showCategory = 10) + 
  ggtitle("Top 10 GO Enrichment Tissue LiverZT12")
dev.off()
#ZT12vsZT12





```

#### 2. KEGG Pathway


```rscript
# Map gene symbols to Entrez IDs
degenes_map <- mapIds(org.Mm.eg.db,
                      keys = deg_genes,
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")


degenes_map <- na.omit(degenes_map)


# KEGG Enrichment Analysis
kegg_results <- enrichKEGG(gene = degenes_map,
                           organism = "mmu",
                           pvalueCutoff = 0.1) 


# Bar plot for KEGG enrichment results
png(filename = "12tissRplot10kegg.png", width = 800, height = 700, res = 90)
dotplot(kegg_results, showCategory = 10) + 
  ggtitle("Top 10 KEGG Pathways issue LiverZT12")
dev.off()
```



```rscript
# Convert results to data frames
go_summary <- as.data.frame(go_results)
kegg_summary <- as.data.frame(kegg_results)

# Display the summaries
head(go_summary)
head(kegg_summary)

write.csv(go_summary, "go_enrichment_tissueLZT12.csv", row.names = FALSE)
write.csv(kegg_summary, "kegg_enrichment_tissueLZT12.csv", row.names = FALSE)
```
