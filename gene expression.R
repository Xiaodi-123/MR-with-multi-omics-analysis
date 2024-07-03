scRNA=readRDS('./scRNA_anno.RDS')
load('table1.Rdata')
gene=unique(table1$Gene)
library(Seurat)
library(viridis)
DotPlot(scRNA,features = gene,cols = c('#dadada','#bc3c29'))
scRNA_T=readRDS('./scRNA_T.RDS')
FeaturePlot(scRNA_T,features = 'N4BP2L1',label = T,pt.size = 0.5,order = T,cols = c('#dadada','#bc3c29'))

