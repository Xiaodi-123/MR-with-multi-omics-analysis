
devtools::install_github("SGDDNB/GeneSwitches")
install.packages('fastglm')
library(GeneSwitches)
library(SingleCellExperiment)
library(mixtools)


scRNA_cm=subset(scRNA_T,T_celltype=='CD8_EFF')

cellinfo <- scRNA_cm@meta.data

sce <- as.SingleCellExperiment(scRNA_cm)

library(slingshot)
library(RColorBrewer)
library(Seurat)

sce_slingshot <- slingshot(sce , clusterLabels = 'T_celltype', reducedDim = 'UMAP', 
                           start.clus = c(3,5), shrink = 0.2)

dev.off()
cl1 <- cellinfo$T_celltype
plot(reducedDims(sce_slingshot)$UMAP,col = brewer.pal(12,"Paired")[cl1],pch=16,asp=1)

igraph::igraph.options(sparsematrices = FALSE)

lines(SlingshotDataSet(sce_slingshot), lwd=2, type = 'lineages', col = 'black')

legend("right",legend = unique(sce$T_celltype),
       col = unique(brewer.pal(12,"Paired")[cl1]),inset=c(3,2,4), pch = 16)

allexpdata <- as.matrix(scRNA_cm@assays$RNA@data);dim(allexpdata)
allcells<-colData(sce_slingshot);dim(allcells)

allcells$slingPseudotime_1

cells <- allcells[!is.na(allcells$slingPseudotime_1),];dim(cells)
expdata <- allexpdata[,rownames(cells)];dim(expdata)

expdata <- expdata[apply(expdata > 0,1,sum) >= 5,];dim(expdata)

rd_UMAP <- Embeddings(object = scRNA_cm, reduction = "umap");dim(rd_UMAP)#原 object = seu3obj.integrated
rd_UMAP <- rd_UMAP[rownames(cells), ];dim(rd_UMAP)
all(rownames(rd_UMAP) == colnames(expdata))

sce <- SingleCellExperiment(assays = List(expdata = expdata))

colData(sce)$Pseudotime <- cells$slingPseudotime_1

reducedDims(sce) <- SimpleList(UMAP = rd_UMAP)
sce_p1 <- sce

h <- hist(assays(sce_p1)$expdata, breaks = 200, plot = FALSE)
{plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
      xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
  abline(v=0.2, col="blue")}

sce_p1 <- binarize_exp(sce_p1, fix_cutoff = TRUE, binarize_cutoff = 0.2)

sce_p1 <- find_switch_logistic_fastglm(sce_p1, downsample = TRUE, show_warning = FALSE, zero_ratio = 0.65, ds_cutoff = 0.65)

table(rowData(sce_p1)$prd_quality)

sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92);dim(sg_allgenes)

sg_gtypes <- filter_switchgenes(sce_p1, allgenes = FALSE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92,
                                genelists = gs_genelists);dim(sg_gtypes)#, genetype = c("Surface proteins", "TFs"))

sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),]);dim(sg_vis)

gl <- unique(table1$Gene)
intersect(sg_vis$geneID, gl)
sg_my <- rowData(sce_p1)[gl,];head(sg_my)
sg_my$feature_type <- "Mendelian genes"
sg_vis <- rbind(sg_vis, sg_my)
plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3.5)

a=sce_p1@assays@data$expdata['N4BP2L1',]
b=sce_p1$Pseudotime

df=data.frame(gene=a,time=b)

ggstatsplot::ggscatterstats(data=df,x='time',y='gene')
