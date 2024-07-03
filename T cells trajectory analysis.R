scRNA_T=readRDS('./scRNA_T.RDS')
BiocManager::install('slingshot')
BiocManager::install('SingleCellExperiment')
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)
table(scRNA_T$T_celltype)
cellinfo <- scRNA_T@meta.data
sce <- as.SingleCellExperiment(scRNA_T)
sce_slingshot <- slingshot(sce , clusterLabels = 'T_celltype', reducedDim = 'UMAP', 
                           start.clus = c(3,5), shrink = 0.2)
cl1 <- cellinfo$T_celltype
plot(reducedDims(sce_slingshot)$UMAP,col = brewer.pal(12,"Paired")[cl1],pch=16,asp=1)
igraph::igraph.options(sparsematrices = FALSE)
lines(SlingshotDataSet(sce_slingshot), lwd=2, type = 'lineages', col = 'black')
legend("right",legend = unique(sce$T_celltype),
       col = unique(brewer.pal(12,"Paired")[cl1]),inset=c(0.5,2,4), pch = 16)
