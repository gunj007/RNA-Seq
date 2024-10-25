# Function Enrichment Analysis
> Enrichment analysis in functional genomics helps identify over-represented biological functions or pathways within a list of genes differentially expressed genes, interpreting biological significance.


```rscript
#install and load the required packages
#BiocManager::install("clusterProfiler"")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
```

#### 1. Gene Ontology (GO)

```rscript
res = read.csv("DESeqResults.csv")
head(res)
```

#### 2. KEGG Pathway







