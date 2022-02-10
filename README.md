# GSEA
Gene set enrichment analysis

ClusterProfiler can be used to assess GSEA of existing gene sets and user-defined gene sets.

An example for user-defined gene sets:

```R
library(clusterProfiler)

setwd('~/OneDrive - Imperial College London/Cebola lab OneDrive folder/Dry-lab/AMP-T2D/LSECs/RNA-seq/GSEA/')
PAvsins=read.table('PA-vs-ins-DEGs-LFC-apeglm-CPM1-ALL.txt')
marker=read.delim('LSEC_marker_genes_ENSMBL.txt')
ERG=read.delim('ERG_responsive_genes_HUVEC_siERG_FDR0.05_log2FC1.5.txt')

PAvsins=PAvsins[order(PAvsins$log2FoldChange,decreasing=TRUE),]
geneList=PAvsins$log2FoldChange
names(geneList)=rownames(PAvsins)

gene.sets=as.data.frame(rbind(cbind('LSEC marker genes',marker[,2]),
                        cbind('ERG responsive genes',ERG[,2])))

results=GSEA(geneList,TERM2GENE=gene.sets)

head(summary(results))
gseaplot(results, "LSEC marker genes")
gseaplot(results, "ERG responsive genes")
```

An example using *piano* to test for directional enrichment, e.g.:

```R
library(piano)
gsc<-loadGSC(rbind(cbind(genes1,'name1'),cbind(genes2,'name2')))

results=runGSA(input,geneSetStat='gsea',signifMethod="geneSampling",adjMethod="fdr",gsc=gsc, nPerm=10000, verbose=TRUE)

GSAsummaryTable(results)
```