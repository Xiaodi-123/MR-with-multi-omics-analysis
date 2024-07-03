scRNA=readRDS('./scRNA_anno.RDS')
scRNA_other=subset(scRNA,celltype != 'T_cells')
rm(scRNA)
gc()
scRNA_other$T_celltype =scRNA_other$celltype
scRNA_chat=merge(scRNA_other,scRNA_T)
rm(scRNA_T)
rm(scRNA_other)
gc()
scRNA_chat_ILD=subset(scRNA_chat,tissue_type=='ILD')
set.seed(123)
a=sample(1:ncol(scRNA_chat_ILD),2000)
scRNA_chat_ILD=scRNA_chat_ILD[,a]

meta =scRNA_chat_ILD@meta.data 
data_input <- as.matrix(scRNA_chat_ILD@assays$RNA@data)
identical(colnames(data_input),rownames(meta))
install.packages('devtools')
library(devtools)
devtools::install_github('sqjin/CellChat')
BiocManager::install('ComplexHeatmap')
library('ComplexHeatmap')
library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "T_celltype")
CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 
dplyr::glimpse(CellChatDB$interaction)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net<- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, sources.use = 'CD8_EFF',
                 title.name = "Number of interactions")
dev.off()
p_bubble= netVisual_bubble(cellchat,
                           sources.use = 'CD8_EFF',
                           remove.isolate = FALSE)+coord_flip()
p_bubble
scRNA_chat_IPF=subset(scRNA_chat,tissue_type=='IPF')
set.seed(123)
a=sample(1:ncol(scRNA_chat_IPF),2000)
scRNA_chat_IPF=scRNA_chat_IPF[,a]
meta =scRNA_chat_IPF@meta.data 
gc()
data_input <- as.matrix(scRNA_chat_IPF@assays$RNA@data)
identical(colnames(data_input),rownames(meta))
library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "T_celltype")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, sources.use = 'CD8_EFF',
                 title.name = "Number of interactions")
dev.off()

pp_bubble= netVisual_bubble(cellchat,
                            sources.use = 'CD8_EFF',
                            remove.isolate = FALSE)+coord_flip()
pp_bubble

