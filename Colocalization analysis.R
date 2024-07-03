load('mr_input_res.Rdata')
library(vcfR)
data <- vcfR::read.vcfR("./eqtl-a-ENSG00000139597/eqtl-a-ENSG00000139597.vcf")
gt=data@gt
gt=as.data.frame(gt)

colnames(gt)
gt$FORMAT[1]

library(tidyverse)

gt$`eqtl-a-ENSG00000239713`[1]
gt=separate(gt,col='eqtl-a-ENSG00000139597',into = c('ES', 'SE',
                                                     'LP','AF','SS',
                                                     'ID'),sep = '\\:')

gc()
gt=na.omit(gt)
colnames(gt)=c('format','beta','se','logpvalue','eaf','samplesize','snp')
gt$beta=as.numeric(gt$beta)
gt$se=as.numeric(gt$se)
gt$logpvalue=as.numeric(gt$logpvalue)
gt$eaf=as.numeric(gt$eaf)
gt$samplesize=as.numeric(gt$samplesize)
gc()
gt$format=NULL
fix=data@fix
fix=as.data.frame(fix)
colnames(fix)
colnames(fix)=c('chr','pos','snp','ref','alt')
fix=fix[,1:5]
eqtl=left_join(fix,gt,by='snp')
eqtl=na.omit(eqtl)
eqtl$maf = ifelse(eqtl$eaf < 0.5, 
                  eqtl$eaf,
                  1 - eqtl$eaf)
eqtl$eaf=NULL
eqtl=eqtl[eqtl$chr==13,]
eqtl$logpvalue=as.numeric(eqtl$logpvalue)
eqtl$p_value=10^(-eqtl$logpvalue)

eqtl$pos=as.numeric(eqtl$pos)
eqtl=eqtl[eqtl$pos > 32997745-1000000 ,]
eqtl=eqtl[eqtl$pos < 32997745+1000000 ,]
my_eqtl=eqtl[,c('snp','p_value','maf')]

colnames(my_eqtl)=c('snp','pvalues','MAF')
my_eqtl=na.omit(my_eqtl)

my_eqtl=my_eqtl[my_eqtl$MAF>0 ,]

library(TwoSampleMR)
coloc_ILD_dat <- extract_outcome_data(snps = c(my_eqtl$snp),
                                      outcomes = "ebi-a-GCST90018643",
                                      proxies = F) %>% 
  mutate(chr.outcome = as.numeric(chr),
         pos.outcome = as.numeric(pos),
         outcome = "Interstitial lung disease", 
         id.outcome = "ebi-a-GCST90018643")

gwas=coloc_ILD_dat
gwas$beta=as.numeric(gwas$beta.outcome)
gwas$se=as.numeric(gwas$se.outcome)
gwas$varbeta=(gwas$se)^2

gwas=gwas[,c('SNP','pval.outcome',"beta",'varbeta')]

colnames(gwas)=c('snp','pvalues','beta','varbeta')

gwas=na.omit(gwas)

library(coloc)
input <- merge(my_eqtl, gwas, by="snp", all=FALSE, suffixes=c("_eqtl","_gwas"))
input=input[!duplicated(input$snp),]
coloc.abf(dataset1=list(pvalues=input$pvalues_gwas, type="cc", s=0.048, N=469827),
          dataset2=list(pvalues=input$pvalues_eqtl, type="quant", N=26279),
          MAF=input$MAF)


