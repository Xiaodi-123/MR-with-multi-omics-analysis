install.packages("Seurat")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("patchwork")
install.packages("hdf5r")
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(hdf5r)

current_directory <- getwd()
subfolder <- "data"
new_path <- file.path(current_directory, subfolder)
setwd(new_path)
rm(current_directory)
rm(new_path)
rm(subfolder)

h5_file <- "./CT1/GSM3489182_Donor_01_raw_gene_bc_matrices_h5.h5"
h5_file2 <- "./CT2/GSM3489185_Donor_02_raw_gene_bc_matrices_h5.h5"
h5_file3<- "./IPF1/GSM3489183_IPF_01_raw_gene_bc_matrices_h5.h5"
h5_file4<- "./IPF2/GSM3489184_IPF_02_raw_gene_bc_matrices_h5.h5"
h5_file5<- "./SSc-ILD1/GSM3489194_SSc-ILD_01_raw_gene_bc_matrices_h5.h5"
h5_file6<- "./SSc-ILD2/GSM3489198_SSc-ILD_02_raw_gene_bc_matrices_h5.h5"

CT1 <- Read10X_h5(file = h5_file)
CT2<- Read10X_h5(file = h5_file2)
IPF1<- Read10X_h5(file = h5_file3)
IPF2<- Read10X_h5(file = h5_file4)
ILD1<- Read10X_h5(file = h5_file5)
ILD2<- Read10X_h5(file = h5_file6)
rm(h5_file)
rm(h5_file2)
rm(h5_file3)
rm(h5_file4)
rm(h5_file5)
rm(h5_file6)
gc()

IPF1=CreateSeuratObject(IPF1,project = 'IPF1')
IPF2=CreateSeuratObject(IPF2,project = 'IPF2')
CT1=CreateSeuratObject(CT1,project = 'CT1')
CT2=CreateSeuratObject(CT2,project = 'CT2')
ILD1=CreateSeuratObject(ILD1,project = 'ILD1')
ILD2=CreateSeuratObject(ILD2,project = 'ILD2')
scRNA <- merge(IPF1,y=c(IPF2,CT1,CT2,ILD1,ILD2))
gc()
rm(CT1)
rm(CT2)
rm(ILD1)
rm(ILD2)
rm(IPF1)
rm(IPF2)

scRNA@meta.data$tissue_type=scRNA@meta.data$orig.ident
scRNA@meta.data$tissue_type=stringr::str_remove(scRNA@meta.data$tissue_type,'[0-9]')

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
library(ggplot2)
col.num <- length(levels(scRNA@active.ident))
VlnPlot(scRNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols =rainbow(col.num), 
        pt.size = 0.01, 
        ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
print(c("Please enter the allowed gene number and ribosome ratio, as shown in the example below：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=200
maxGene=4000
pctMT=10
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(scRNA@active.ident))
VlnPlot(scRNA,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols =rainbow(col.num), 
        pt.size = 0.1, 
        ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")
library(harmony)
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")
DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNA_harmony,reduction = 'harmony')
scRNA <- FindNeighbors(scRNA_harmony, dims = 1:10, reduction = "harmony")
scRNA <- FindClusters(scRNA)
scRNA <- RunUMAP(scRNA, dims = 1:10,reduction = 'harmony')
DimPlot(scRNA,split.by = 'tissue_type')
rm(scRNA_harmony)
install.packages("BiocManager",force=TRUE)
BiocManager::install("SingleR",force = TRUE)
library(SingleR)
BiocManager::install("celldex",force = TRUE)
install.packages("RCurl")
library(RCurl)
install.packages("GenomeInfoDb")
library(celldex)
refdata <- celldex::HumanPrimaryCellAtlasData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata,
                    labels =refdata$label.main,
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
Idents(scRNA)=scRNA$celltype
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
DimPlot(scRNA,split.by = 'tissue_type',cols = col_vector[20:30])
saveRDS(scRNA,file = 'scRNA_anno.RDS')
scRNA_T=subset(scRNA,celltype=='T_cells')
scRNAsub<- FindVariableFeatures(scRNA_T, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNAsub)
library(harmony)
set.seed(1000)
scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")
ElbowPlot(scRNAsub,reduction = 'harmony')
scRNAsub <- FindNeighbors(scRNAsub, reduction = 'harmony',dims = 1:10)
scRNAsub <- FindClusters(scRNAsub)
scRNAsub <- RunUMAP(scRNAsub, reduction = 'harmony',dims = 1:10)
DimPlot(scRNAsub, label=TRUE,split.by = 'tissue_type') 
dev.off()
library(reshape2)
library(ggplot2)
prop_df <- table(scRNAsub@meta.data$seurat_clusters,scRNAsub@meta.data$tissue_type) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
sample_color <- col_vector[1:10] 
prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
prop
Idents(scRNAsub)=scRNAsub$seurat_clusters
rm(T_marker)
T_marker <- c('KLRB1','IL7R','CXCR6','CD4',#CD4_RM
              'CXCR5','BCL6','ICA1','CD200','PDCD1', #CD4_FH
              'IL23R',"CCR6",'CAPG','RORC','IL17A', #TH17
              'FOXP3','IL2RA','IL1R2',#CD4_REG
              'CXCL13','IFNG','CXCR3','ICOS','GZMB',#TH1
              'GATA3','IL13','IL1RL1',#TH2
              'CD8A','GZMK',  # CD8_EM
              'GZMA','CCL5',  #CD8_CM
              'CD6','XCL1','XCL2','CD69','RORA', #CD8_RM
              'GZMH','S1PR5','KLRG1','CX3CR1','PRF1',#CD8_EFF
              'HAVCR2','LAG3', # CD8_exhau
              'EPCAM','CD19','CD3E') 
DotPlot(scRNAsub, 
        features = T_marker,
        group.by = "seurat_clusters") + coord_flip()
T_celltype=c('CD4_RM',
             'T cell',
             'T cell',
             'CD8_CM',
             'CD8_EM',
             'CD8_EFF',
             'CD4_REG',
             'NK',
             'CD8_RM')
Idents(scRNAsub) <- scRNAsub@meta.data$seurat_clusters
names(T_celltype) <- levels(scRNAsub)
scRNAsub<- RenameIdents(scRNAsub, T_celltype)

scRNAsub@meta.data$T_celltype <- Idents(scRNAsub)
Idents(scRNAsub)=scRNAsub@meta.data$T_celltype

colors=c('#313c63','#b42e20','#ebc03e','#377b4c',
         '#7bc7cd','#5d84a4','#bc3c29')
DimPlot(scRNAsub, group.by="T_celltype", label=T, label.size=5,cols = col_vector[30:50],
        pt.size = 1,split.by = 'tissue_type')
library(reshape2)
library(ggplot2)
prop_df <- table(scRNAsub@meta.data$T_celltype,scRNAsub@meta.data$tissue_type) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

sample_color <- col_vector[1:10] 

prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
prop
saveRDS(scRNAsub,file = 'scRNA_T.RDS')
