library(data.table)
rt=fread('./GSE231693_count1.tsv',data.table = F)
rownames(rt)=rt$V1
rt$V1=NULL
devtools::install_github("IOBR/IOBR")
BiocManager::install('preprocessCore')
BiocManager::install('biomaRt')
BiocManager::install('DESeq2')
BiocManager::install('limma')
library(IOBR)
rm(scRNA_chat_ILD)
rm(scRNA_EFF)
rm(scRNA_T)
gc()
rt=count2tpm(rt,idType = 'Ensembl',org = 'hsa')

rt=as.data.frame(rt)

max(rt)
rt=log2(rt+1)

a1=grep('SSc',colnames(rt))

exp1=rt[,a1]
exp2=rt[,-a1]

rt=cbind(exp2,exp1)


load('table1.Rdata')

data=rt[unique(table1$Gene),]

anno=data.frame(row.names =colnames(rt),group=c(rep('Healthy',20),
                                                rep('SSc',20)))
pheatmap::pheatmap(data,cluster_cols = F,
                   scale = 'row',show_colnames = F,annotation_col = anno)

df=data.frame(gene=as.numeric(rt['N4BP2L1',]),group=anno$group)
ggpubr::ggboxplot(data = df, x = 'group',y='gene',color = 'group',palette = 'jco',notch = T,size = 1)+
  stat_compare_means()

