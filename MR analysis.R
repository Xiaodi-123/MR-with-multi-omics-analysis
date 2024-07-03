library(dplyr)
scRNA=readRDS('./scRNA_anno.RDS')
scRNAsub=readRDS('./scRNA_T.RDS')
library(Seurat)
scRNA_other=subset(scRNA,celltype !='T_cells')
scRNA_EFF=subset(scRNAsub,T_celltype == 'CD8_EFF')
gc()
scRNA_other$T_celltype=scRNA_other$celltype
scRNA_compare=merge(scRNA_other,scRNA_EFF)
table(scRNA_compare$T_celltype)
gc()
rm()

df_EFF=FindMarkers(scRNAsub,ident.1 = 'CD8_EFF',only.pos = T,logfc.threshold = 0.5)

df_T=FindMarkers(scRNA_compare,ident.1 = 'CD8_EFF',only.pos = T,logfc.threshold = 0.5)

ss=intersect(rownames(df_EFF),rownames(df_T))

ss

save(ss,file ='key_marker_gene.Rdata')

load('key_marker_gene.Rdata')
write.csv(ss,file ='key_marker_gene.csv',quote = F)

gene=ss
gene=gene[!duplicated(gene)]

BiocManager::install('clusterProfiler')
library(clusterProfiler)
BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
df=bitr(geneID = gene,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
id=df$ENSEMBL
id=id[!duplicated(id)]
id

devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
https://gwas.mrcieu.ac.uk/datasets/
  
  exposure_id=paste0('eqtl-a-',id)
exposure_dat <- extract_instruments(exposure_id, p1=5e-08, clump=TRUE)

R2a=2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2b=2*exposure_dat$se.exposure*exposure_dat$se.exposure*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
R2=R2a/(R2a+R2b)

exposure_dat$F_statistics=R2*(exposure_dat$samplesize.exposure-2)/(1-R2)

exposure_dat=exposure_dat[exposure_dat$F_statistics>10,]

exposure_dat$exposure=stringr::str_sub(exposure_dat$exposure,1,15)

exposure_dat$ENSEMBL =exposure_dat$exposure

library(tidyverse)
exposure_dat= left_join(exposure_dat,df, by = "ENSEMBL")

exposure_dat$exposure=NULL
exposure_dat$exposure=exposure_dat$SYMBOL
exposure_dat$SYMBOL=NULL
save(exposure_dat,file ='exposure_dat.Rdata')

load('exposure_dat.Rdata')

rm(scRNA_CM)
rm(scRNA_compare)
rm(scRNA)
gc()

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST90018643")
exposure_dat=exposure_dat[exposure_dat$SNP %in% outcome_dat$SNP,]
harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)

mr_modified <- function(dat = harmonised_dat, prop_var_explained = T)
{
  mr_res <- mr(dat)
  
  pve <- dat %>% 
    dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  
  return(mr_res)
}

mr_res <- mr_modified(harmonised_dat, prop_var_explained = T)

save(exposure_dat,outcome_dat,mr_res,harmonised_dat,file ='mr_input_res.Rdata')

load('mr_input_res.Rdata')

heter_dat=mr_heterogeneity(harmonised_dat)
write.csv(heter_dat, file="table.heterogeneity.csv", row.names=F)

pleio_dat=mr_pleiotropy_test(harmonised_dat)
write.csv(pleio_dat, file="table.pleiotropy.csv", row.names=F)

table1 <- mr_res %>% 
  filter(pval < 0.05,
         method %in% c("Wald ratio","Inverse variance weighted")) %>% 
  left_join(exposure_dat, by = "exposure")

table1 <- table1 %>% 
  generate_odds_ratios()%>% 
  mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",or,or_lci95, or_uci95),
         `P value` = scales::scientific(pval),
         `PVE` = paste0(sprintf("%.2f",100*pve),"%"),
         `F statistics` = sprintf("%.2f",F_statistics)) %>% 
  dplyr::select(Gene = exposure, `ENSEMBL ID` = ENSEMBL,
                SNP, `Effect allele` = effect_allele.exposure, 
                `OR (95% CI)`, `P value`, 
                PVE, `F statistics`)

save(table1,file ='table1.Rdata')

volcano_plot <- function(.data, 
                         number_comparasion = 1,
                         title = "eQTL",
                         col_beta = "b",
                         col_size = "pve",
                         col_label = "exposure",
                         legend.position = "none")
{
  p_thershold <- 0.05/number_comparasion
  
  p <- .data %>% 
    rename(beta := !!col_beta,
           size := !!col_size,
           label := !!col_label) %>% 
    mutate(x = beta,
           y = -log10(pval),
           label = ifelse(pval < p_thershold, label, NA)) %>% 
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(size = size), alpha = 0.5, color = "#0072b5") +
    geom_vline(xintercept = 0, linetype = 2)+
    geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6.5),
          legend.position = legend.position)+
    labs(x = "ln(OR)", 
         y = parse(text = "-log[10]*(italic(P)-value)"),
         title = title) +
    scale_size(name = "PVE",
               breaks = c(0.2*1:3)) +
    ggrepel::geom_label_repel(aes(label = label),size = 3)
  plot(p)
}



mr_res %>% 
  filter(method %in% c("Wald ratio","Inverse variance weighted")) %>% 
  volcano_plot(number_comparasion = 1)

mr_res1=mr_res[mr_res$exposure=='GZMH',]
harmonised_dat1=harmonised_dat[harmonised_dat$exposure=='GZMH',]

mr_scatter_plot(mr_res1, harmonised_dat1)

res_single=mr_singlesnp(harmonised_dat1)
mr_forest_plot(res_single)

mr_funnel_plot(singlesnp_results = res_single)

mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonised_dat1))


mr_res2=mr_res[mr_res$exposure %in% table1$Gene,]
result_or=generate_odds_ratios(mr_res2)
write.table(result_or[,4:ncol(result_or)],"OR.txt",row.names = F,sep = "\t",quote = F)

library(grid)
library(forestploter)
mydata=read.table("OR.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
forest(mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)