scRNA_T=readRDS('./scRNA_T.RDS')

gc()
scRNA_EFF=subset(scRNA_T,T_celltype=='CD8_EFF')
scRNA_EFF$gene_group=ifelse(scRNA_EFF@assays$RNA@counts['N4BP2L1',]>0,'N4BP2L1+CD8_EFF','N4BP2L1-CD8_EFF')

scRNA_otherT=subset(scRNA_T,T_celltype != 'CD8_EFF')

# 加列
scRNA_otherT$gene_group=scRNA_otherT$T_celltype

scRNA_metab=merge(scRNA_EFF,c(scRNA_otherT))

gc()

rm(scRNA_chat)
gc()
scRNA_metab_ILD=subset(scRNA_metab,tissue_type=='ILD')

set.seed(123)
a=sample(1:ncol(scRNA_metab_ILD),2000)
scRNA_metab_ILD=scRNA_metab_ILD[,a]

devtools::install_github("YosefLab/VISION")
BiocManager::install('AUCell')
BiocManager::install('GSVA')
library(AUCell)
devtools::install_github("wu-yc/scMetabolism")
library(scMetabolism)
library(ggplot2)
library(rsvd)
scRNA_metab_ILD<-sc.metabolism.Seurat(obj = scRNA_metab_ILD, method = 'AUCell', imputation = F, ncores = 2, metabolism.type = "KEGG")
input.pathway <- rownames(scRNA_metab_ILD@assays[["METABOLISM"]][["score"]])[61:90]
DotPlot.metabolism(obj =scRNA_metab_ILD,
                   pathway = input.pathway, phenotype = "gene_group", norm = "y")
