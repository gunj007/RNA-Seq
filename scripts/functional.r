library(enrichplot)
library(ggplot2)
library(clusterProfiler)
library(pathview)


restiss <- results(ds, name ="typetissue_liver_vs_heart")
resint <- results(ds, name ="typetissueliver.timeZT12")
resime <- results(ds, name = "time_ZT12_vs_ZT0")

require(DOSE)
graphics.off()
# rm na value
deg_genes <- rownames(restiss[!is.na(restiss$padj) & restiss$padj < 0.05, ])
deg_genes <- rownames(restime[!is.na(restime$padj) & restime$padj < 0.05, ])
deg_genes <- rownames(resint[!is.na(resint$padj) & resint$padj < 0.05, ])
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
png(filename = "12tissRplot10go.png", width = 800, height = 700, res = 100)
dotplot(go_results, showCategory = 10) + 
  ggtitle("Top 10 GO Enrichment Tissue Specific")
dev.off()
#ZT12vsZT12

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
png(filename = "13tissRplot10kegg.png", width = 800, height = 700, res = 90)
dotplot(kegg_results, showCategory = 10) + 
  ggtitle("Top 10 KEGG Pathways Tissue Specific")
dev.off()

# Convert results to data frames
go_summary <- as.data.frame(go_results)
kegg_summary <- as.data.frame(kegg_results)

# Display the summaries
head(go_summary)
head(kegg_summary)

write.csv(go_summary, "go_enrichment_tiss.csv", row.names = FALSE)
write.csv(kegg_summary, "kegg_enrichment_tiss.csv", row.names = FALSE)




